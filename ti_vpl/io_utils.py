from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def element_centroids(points_m: np.ndarray, tets: np.ndarray) -> np.ndarray:
    return points_m[tets].mean(axis=1)


def spherical_roi(centroids_m: np.ndarray, center_mm: np.ndarray, radius_mm: float) -> np.ndarray:
    center_m = np.asarray(center_mm, dtype=float) * 1e-3
    return np.linalg.norm(centroids_m - center_m[None, :], axis=1) <= float(radius_mm) * 1e-3


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values)
    v = np.asarray(values)[order]
    w = np.asarray(weights)[order]
    cumulative = np.cumsum(w) / np.maximum(np.sum(w), 1e-30)
    return float(v[np.searchsorted(cumulative, quantile, side="left")])


def roi_statistics(values: np.ndarray, volumes_m3: np.ndarray, mask: np.ndarray) -> dict[str, float | int]:
    v = np.asarray(values)[mask]
    w = np.asarray(volumes_m3)[mask]
    if not len(v):
        raise ValueError("ROI contains no tetrahedra")
    return {
        "tetrahedra": int(len(v)),
        "volume_mm3": float(w.sum() * 1e9),
        "mean": float(np.average(v, weights=w)),
        "median": weighted_quantile(v, w, 0.50),
        "p05": weighted_quantile(v, w, 0.05),
        "p95": weighted_quantile(v, w, 0.95),
        "maximum": float(v.max()),
    }


def write_json(path: str | Path, payload: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

