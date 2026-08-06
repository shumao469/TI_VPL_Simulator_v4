"""Coordinate conventions shared by the web preview and FEM scripts.

The default stereotaxic convention is:
  AP: anterior is positive, posterior is negative
  ML: right is positive, left is negative
  DV: dorsal is positive, ventral/deep is negative

All values in this module are millimetres relative to bregma.  Allen CCF voxel
coordinates are a different coordinate system and must not be substituted
without an explicit registration transform.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class BregmaPoint:
    AP: float
    ML: float
    DV: float

    def as_array(self) -> np.ndarray:
        return np.asarray([self.AP, self.ML, self.DV], dtype=float)

    def as_dict(self) -> dict[str, float]:
        return {k: float(v) for k, v in asdict(self).items()}

    @classmethod
    def from_iterable(cls, values: Iterable[float]) -> "BregmaPoint":
        ap, ml, dv = [float(v) for v in values]
        return cls(ap, ml, dv)


DEFAULT_VPL = BregmaPoint(AP=-1.94, ML=-1.75, DV=-3.10)


def crossed_montage(
    target: BregmaPoint = DEFAULT_VPL,
    half_spacing_mm: float = 2.08,
    electrode_dv_mm: float = 0.50,
    rotation_deg: float = 0.0,
) -> dict[str, tuple[BregmaPoint, BregmaPoint]]:
    """Return the Xu Lab target-centred crossed four-electrode montage.

    At zero rotation, channel 1 is E1--E2 along AP and channel 2 is E3--E4
    along ML.  With the default target/spacing the centres are exactly
    E1=(+0.14,-1.75,+0.50), E2=(-4.02,-1.75,+0.50),
    E3=(-1.94,+0.33,+0.50), E4=(-1.94,-3.83,+0.50) mm from bregma.
    These are initial requested surface coordinates; formal FEM projects them
    to a registered head surface and reports the projection residual.
    """
    theta = np.deg2rad(float(rotation_deg))
    ap_axis = np.asarray([np.cos(theta), np.sin(theta)])
    ml_axis = np.asarray([-np.sin(theta), np.cos(theta)])
    center = np.asarray([target.AP, target.ML])

    def point(offset: np.ndarray) -> BregmaPoint:
        ap, ml = center + float(half_spacing_mm) * offset
        return BregmaPoint(float(ap), float(ml), float(electrode_dv_mm))

    e1, e2 = point(ap_axis), point(-ap_axis)
    e3, e4 = point(ml_axis), point(-ml_axis)
    return {"channel_1": (e1, e2), "channel_2": (e3, e4)}


def apply_affine_mm(points_mm: np.ndarray, affine_4x4: np.ndarray) -> np.ndarray:
    """Apply an explicit millimetre-to-millimetre homogeneous transform."""
    points = np.atleast_2d(np.asarray(points_mm, dtype=float))
    affine = np.asarray(affine_4x4, dtype=float)
    if affine.shape != (4, 4):
        raise ValueError("affine_4x4 must have shape (4, 4)")
    hom = np.c_[points, np.ones(len(points))]
    return (affine @ hom.T).T[:, :3]
