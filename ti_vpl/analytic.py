"""Fast homogeneous-volume preview used by the Streamlit app.

This module deliberately labels its output as a preview.  It is useful for
checking coordinates, montage geometry, TI formulae and interface behaviour,
but it is not a replacement for a segmented, registered FEM head model.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .coordinates import BregmaPoint, DEFAULT_VPL, crossed_montage
from .envelope import maximum_orientation_envelope


@dataclass
class PreviewParameters:
    target: BregmaPoint = DEFAULT_VPL
    channel_1_current_mA: float = 0.50
    channel_2_current_mA: float = 0.50
    f1_Hz: float = 2000.0
    f2_Hz: float = 2010.0
    roi_radius_mm: float = 0.60
    half_spacing_mm: float = 2.08
    electrode_dv_mm: float = 0.50
    rotation_deg: float = 0.0
    conductivity_S_per_m: float = 0.275
    skull_transfer_factor: float = 0.028


def is_brain(points_mm: np.ndarray) -> np.ndarray:
    p = np.asarray(points_mm, dtype=float)
    return ((p[..., 0] + 1.25) / 4.95) ** 2 + (p[..., 1] / 4.8) ** 2 + ((p[..., 2] + 3.0) / 3.05) ** 2 <= 1.0


def point_pair_field(
    points_mm: np.ndarray,
    source_mm: np.ndarray,
    sink_mm: np.ndarray,
    current_mA: float,
    conductivity_S_per_m: float = 0.275,
    skull_transfer_factor: float = 0.028,
) -> np.ndarray:
    """Analytic field from a source/sink pair in a homogeneous conductor."""
    p = np.asarray(points_mm, dtype=float) * 1e-3
    source = np.asarray(source_mm, dtype=float) * 1e-3
    sink = np.asarray(sink_mm, dtype=float) * 1e-3
    current_A = float(current_mA) * 1e-3
    gain = current_A * float(skull_transfer_factor) / (4.0 * np.pi * float(conductivity_S_per_m))

    def contribution(electrode: np.ndarray, sign: float) -> np.ndarray:
        rvec = p - electrode
        r = np.maximum(np.linalg.norm(rvec, axis=-1), 0.45e-3)
        return sign * gain * rvec / r[..., None] ** 3

    return contribution(source, +1.0) + contribution(sink, -1.0)


def preview_fields(points_mm: np.ndarray, params: PreviewParameters) -> dict[str, np.ndarray]:
    m = crossed_montage(
        params.target,
        params.half_spacing_mm,
        params.electrode_dv_mm,
        params.rotation_deg,
    )
    ch1 = np.vstack([p.as_array() for p in m["channel_1"]])
    ch2 = np.vstack([p.as_array() for p in m["channel_2"]])
    common = dict(
        conductivity_S_per_m=params.conductivity_S_per_m,
        skull_transfer_factor=params.skull_transfer_factor,
    )
    e1 = point_pair_field(points_mm, ch1[0], ch1[1], current_mA=params.channel_1_current_mA, **common)
    e2 = point_pair_field(points_mm, ch2[0], ch2[1], current_mA=params.channel_2_current_mA, **common)
    modulation = maximum_orientation_envelope(e1, e2)
    if np.isclose(params.f1_Hz, params.f2_Hz):
        modulation = np.zeros_like(modulation)
    return {"E1": e1, "E2": e2, "mTI": modulation}


def orthogonal_slice(kind: str, params: PreviewParameters, n: int = 90) -> dict[str, np.ndarray]:
    if kind not in {"coronal", "sagittal", "axial"}:
        raise ValueError("kind must be coronal, sagittal, or axial")
    target = params.target
    if kind == "coronal":
        x = np.linspace(-5.2, 5.2, n)   # ML
        y = np.linspace(-6.1, 0.15, n) # DV
        xx, yy = np.meshgrid(x, y)
        points = np.stack([np.full_like(xx, target.AP), xx, yy], axis=-1)
        xlabel, ylabel = "ML (mm)", "DV (mm)"
    elif kind == "sagittal":
        x = np.linspace(-6.2, 3.8, n)   # AP
        y = np.linspace(-6.1, 0.15, n) # DV
        xx, yy = np.meshgrid(x, y)
        points = np.stack([xx, np.full_like(xx, target.ML), yy], axis=-1)
        xlabel, ylabel = "AP (mm)", "DV (mm)"
    else:
        x = np.linspace(-5.2, 5.2, n)  # ML
        y = np.linspace(-6.2, 3.8, n)  # AP
        xx, yy = np.meshgrid(x, y)
        points = np.stack([yy, xx, np.full_like(xx, target.DV)], axis=-1)
        xlabel, ylabel = "ML (mm)", "AP (mm)"
    values = preview_fields(points, params)["mTI"]
    values = np.where(is_brain(points), values, np.nan)
    return {"x": x, "y": y, "z": values, "points": points, "xlabel": xlabel, "ylabel": ylabel}


def preview_metrics(params: PreviewParameters, shape: tuple[int, int, int] = (34, 35, 28)) -> dict[str, float]:
    ap = np.linspace(-6.0, 3.5, shape[0])
    ml = np.linspace(-5.0, 5.0, shape[1])
    dv = np.linspace(-6.0, 0.0, shape[2])
    aa, mm, dd = np.meshgrid(ap, ml, dv, indexing="ij")
    points = np.stack([aa, mm, dd], axis=-1)
    brain = is_brain(points)
    fields = preview_fields(points, params)["mTI"]
    target = params.target.as_array()
    roi = np.linalg.norm(points - target, axis=-1) <= params.roi_radius_mm
    roi_values = fields[brain & roi]
    brain_values = fields[brain]
    if not len(roi_values):
        raise ValueError("The requested ROI contains no preview grid samples")

    threshold = np.quantile(brain_values, 0.75)
    high = brain & (fields >= threshold)
    weights = np.where(high, fields, 0.0)
    centroid = np.sum(points * weights[..., None], axis=(0, 1, 2)) / np.maximum(weights.sum(), 1e-30)
    return {
        "vpl_mean_V_per_m": float(np.mean(roi_values)),
        "vpl_median_V_per_m": float(np.median(roi_values)),
        "vpl_p95_V_per_m": float(np.quantile(roi_values, 0.95)),
        "brain_p95_V_per_m": float(np.quantile(brain_values, 0.95)),
        "brain_p99_V_per_m": float(np.quantile(brain_values, 0.99)),
        "focality_median_over_brain_p95": float(np.median(roi_values) / np.maximum(np.quantile(brain_values, 0.95), 1e-30)),
        "high_field_centroid_distance_mm": float(np.linalg.norm(centroid - target)),
    }
