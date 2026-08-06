#!/usr/bin/env python
"""Attach an already registered subject-space VPL NIfTI mask to a VTU mesh."""
from __future__ import annotations

import argparse
from pathlib import Path

import meshio
import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", default="outputs/ti_vpl_fields.vtu")
    parser.add_argument("--vpl-mask", required=True)
    parser.add_argument("--output", default="outputs/ti_vpl_fields_registered_roi.vtu")
    parser.add_argument(
        "--mesh-to-nifti-affine",
        help="Optional text file containing a 4x4 transform from mesh mm to NIfTI RAS mm",
    )
    args = parser.parse_args()

    mesh = meshio.read(args.mesh)
    tetra_index = next((i for i, block in enumerate(mesh.cells) if block.type.startswith("tetra")), None)
    if tetra_index is None:
        raise SystemExit("No tetrahedral cells found")
    tets = mesh.cells[tetra_index].data[:, :4]
    centroids_mesh_mm = mesh.points[tets].mean(axis=1)

    transform = np.eye(4) if args.mesh_to_nifti_affine is None else np.loadtxt(args.mesh_to_nifti_affine)
    if transform.shape != (4, 4):
        raise SystemExit("--mesh-to-nifti-affine must contain a 4x4 matrix")
    hom = np.c_[centroids_mesh_mm, np.ones(len(centroids_mesh_mm))]
    centroids_ras_mm = (transform @ hom.T).T[:, :3]

    image = nib.load(args.vpl_mask)
    data = np.asarray(image.dataobj)
    voxel = (np.linalg.inv(image.affine) @ np.c_[centroids_ras_mm, np.ones(len(centroids_ras_mm))].T).T[:, :3]
    sampled = map_coordinates(data.astype(float), voxel.T, order=0, mode="constant", cval=0.0) > 0.5

    cell_data = {name: [np.asarray(values[i]) for i in range(len(mesh.cells))] for name, values in mesh.cell_data.items()}
    cell_data["VPL_registered"] = [
        sampled.astype(np.uint8) if i == tetra_index else np.zeros(len(block.data), dtype=np.uint8)
        for i, block in enumerate(mesh.cells)
    ]
    out = meshio.Mesh(mesh.points, mesh.cells, point_data=mesh.point_data, cell_data=cell_data, field_data=mesh.field_data)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    meshio.write(args.output, out)
    print(f"Wrote {args.output}; VPL tetrahedra={int(sampled.sum())}")


if __name__ == "__main__":
    main()

