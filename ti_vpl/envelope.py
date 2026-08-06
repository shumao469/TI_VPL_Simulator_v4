"""Temporal-interference envelope metrics for two real-valued vector fields."""
from __future__ import annotations

import numpy as np


def maximum_orientation_envelope(e1: np.ndarray, e2: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    """Maximum modulation amplitude over all unit directions.

    Computes max_n 2*min(|n dot E1|, |n dot E2|) independently at each point.
    This is an optimistic orientation-independent upper-bound metric.
    """
    a = np.asarray(e1, dtype=float)
    b = np.asarray(e2, dtype=float)
    if a.shape != b.shape or a.shape[-1] != 3:
        raise ValueError("e1 and e2 must have identical shape (..., 3)")

    na = np.linalg.norm(a, axis=-1)
    nb = np.linalg.norm(b, axis=-1)
    swap = na > nb
    small = np.where(swap[..., None], b, a)
    large = np.where(swap[..., None], a, b)
    A = np.minimum(na, nb)
    B = np.maximum(na, nb)

    dot_abs = np.abs(np.sum(small * large, axis=-1))
    c = np.clip(dot_abs / np.maximum(A * B, eps), 0.0, 1.0)
    aligned = dot_abs >= A * A
    out = np.zeros_like(A)
    out[aligned] = 2.0 * A[aligned]

    other = (~aligned) & (A > eps) & (B > eps)
    cross = A[other] * B[other] * np.sqrt(np.maximum(1.0 - c[other] ** 2, 0.0))
    denom = np.sqrt(np.maximum(A[other] ** 2 + B[other] ** 2 - 2.0 * A[other] * B[other] * c[other], eps))
    out[other] = 2.0 * cross / denom
    return out


def directional_envelope(e1: np.ndarray, e2: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """Modulation amplitude along one direction or one direction per point."""
    n = np.asarray(direction, dtype=float)
    n = n / np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-30)
    p1 = np.abs(np.sum(np.asarray(e1) * n, axis=-1))
    p2 = np.abs(np.sum(np.asarray(e2) * n, axis=-1))
    return 2.0 * np.minimum(p1, p2)


def carrier_rms(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    """Long-time RMS magnitude for two distinct carrier frequencies."""
    return np.sqrt(0.5 * (np.sum(np.abs(e1) ** 2, axis=-1) + np.sum(np.abs(e2) ** 2, axis=-1)))


def vector_angle_deg(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    n1 = np.linalg.norm(e1, axis=-1)
    n2 = np.linalg.norm(e2, axis=-1)
    c = np.sum(e1 * e2, axis=-1) / np.maximum(n1 * n2, 1e-30)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))

