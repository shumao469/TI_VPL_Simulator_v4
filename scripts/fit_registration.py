#!/usr/bin/env python3
"""Fit and QC an affine from bregma AP/ML/DV landmarks to FEM x/y/z."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ti_vpl.registration import fit_affine, write_affine


COLUMNS = ("bregma_AP", "bregma_ML", "bregma_DV", "fem_x", "fem_y", "fem_z")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--landmarks", required=True, help="CSV with bregma_AP,...,fem_z columns")
    parser.add_argument("--output-affine", default="registration/bregma_to_fem_affine.txt")
    parser.add_argument("--output-qc", default="registration/registration_qc.json")
    parser.add_argument("--maximum-rms-mm", type=float, default=0.20)
    args = parser.parse_args()

    with open(args.landmarks, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    missing = [name for name in COLUMNS if not rows or name not in rows[0]]
    if missing:
        raise ValueError(f"Missing CSV columns: {missing}; expected {COLUMNS}")
    values = np.asarray([[float(row[name]) for name in COLUMNS] for row in rows])
    qc = fit_affine(values[:, :3], values[:, 3:])
    payload = {
        "landmark_count": len(values),
        "rms_mm": qc.rms_mm,
        "maximum_mm": qc.maximum_mm,
        "condition_number": qc.condition_number,
        "maximum_allowed_rms_mm": args.maximum_rms_mm,
        "passed": qc.rms_mm <= args.maximum_rms_mm,
        "residual_mm": qc.residual_mm.tolist(),
        "coordinate_order": "input=[AP,ML,DV] bregma mm; output=[x,y,z] FEM mm",
    }
    write_affine(args.output_affine, qc.affine)
    Path(args.output_qc).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_qc).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if not payload["passed"]:
        raise SystemExit("Registration QC failed; improve landmarks before formal FEM")


if __name__ == "__main__":
    main()
