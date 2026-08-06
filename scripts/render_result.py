#!/usr/bin/env python
"""Render a publication-oriented PNG from the FEM VTU file."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="outputs/ti_vpl_fields.vtu")
    parser.add_argument("--scalar", default="TI_envelope_max_V_per_m")
    parser.add_argument("--roi", default="VPL_sphere")
    parser.add_argument("--output", default="outputs/TI_VPL_3D.png")
    parser.add_argument("--percentile", type=float, default=99.0)
    args = parser.parse_args()

    try:
        import pyvista as pv
    except ModuleNotFoundError as exc:
        raise SystemExit("pyvista 未安装：python -m pip install pyvista") from exc

    grid = pv.read(args.input)
    if args.scalar not in grid.cell_data:
        raise SystemExit(f"Scalar not found: {args.scalar}")
    brain = grid.threshold([0.5, 1.5], scalars="brain") if "brain" in grid.cell_data else grid
    values = np.asarray(brain.cell_data[args.scalar])
    clim = [0.0, float(np.quantile(values[np.isfinite(values)], args.percentile / 100.0))]
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 1400))
    plotter.set_background("white")
    plotter.add_mesh(brain, scalars=args.scalar, cmap="turbo", clim=clim, opacity=0.50, show_edges=False, scalar_bar_args={"title": "TI envelope (V/m)"})
    if args.roi in grid.cell_data:
        roi = grid.threshold([0.5, 1.5], scalars=args.roi)
        if roi.n_cells:
            plotter.add_mesh(roi.extract_surface(), color="#fff2a8", opacity=0.95, line_width=3)
    plotter.add_axes(line_width=3, labels_off=False)
    plotter.camera_position = "iso"
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    plotter.show(screenshot=str(output), auto_close=True)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

