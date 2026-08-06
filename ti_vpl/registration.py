"""Bregma stereotaxic to FEM-space registration with explicit QC."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class RegistrationQC:
    affine: np.ndarray
    predicted_mm: np.ndarray
    residual_mm: np.ndarray
    rms_mm: float
    maximum_mm: float
    condition_number: float


def read_affine(path: str | Path) -> np.ndarray:
    matrix = np.loadtxt(path, dtype=float)
    if matrix.shape != (4, 4):
        raise ValueError(f"Registration affine must be 4x4, got {matrix.shape}")
    if not np.all(np.isfinite(matrix)) or abs(np.linalg.det(matrix[:3, :3])) < 1e-10:
        raise ValueError("Registration affine is non-finite or singular")
    if not np.allclose(matrix[3], [0.0, 0.0, 0.0, 1.0], atol=1e-8):
        raise ValueError("Registration affine last row must be [0, 0, 0, 1]")
    return matrix


def write_affine(path: str | Path, affine: np.ndarray) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.asarray(affine, dtype=float), fmt="%.12g")


def apply_affine(points_mm: np.ndarray, affine: np.ndarray) -> np.ndarray:
    points = np.asarray(points_mm, dtype=float)
    scalar = points.ndim == 1
    points = np.atleast_2d(points)
    if points.shape[1] != 3:
        raise ValueError("Points must have shape (N, 3) in AP/ML/DV or FEM x/y/z order")
    homogeneous = np.column_stack([points, np.ones(len(points))])
    transformed = (np.asarray(affine, dtype=float) @ homogeneous.T).T[:, :3]
    return transformed[0] if scalar else transformed


def fit_affine(bregma_ap_ml_dv_mm: np.ndarray, fem_xyz_mm: np.ndarray) -> RegistrationQC:
    source = np.asarray(bregma_ap_ml_dv_mm, dtype=float)
    target = np.asarray(fem_xyz_mm, dtype=float)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3:
        raise ValueError("Source and target landmarks must both have shape (N, 3)")
    if len(source) < 4:
        raise ValueError("At least four non-coplanar landmarks are required for an affine")
    design = np.column_stack([source, np.ones(len(source))])
    rank = np.linalg.matrix_rank(design)
    if rank < 4:
        raise ValueError("Landmarks are coplanar or otherwise rank-deficient")
    coefficients, _, _, _ = np.linalg.lstsq(design, target, rcond=None)
    affine = np.eye(4)
    affine[:3, :] = coefficients.T
    predicted = apply_affine(source, affine)
    residual = np.linalg.norm(predicted - target, axis=1)
    return RegistrationQC(
        affine=affine,
        predicted_mm=predicted,
        residual_mm=residual,
        rms_mm=float(np.sqrt(np.mean(residual**2))),
        maximum_mm=float(np.max(residual)),
        condition_number=float(np.linalg.cond(design)),
    )


def validate_registration(qc: dict, maximum_rms_mm: float) -> None:
    rms = float(qc.get("rms_mm", np.inf))
    if not np.isfinite(rms) or rms > maximum_rms_mm:
        raise RuntimeError(
            f"Registration QC failed: RMS={rms:.3f} mm exceeds {maximum_rms_mm:.3f} mm"
        )
