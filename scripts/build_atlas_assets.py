#!/usr/bin/env python
"""Build compact Allen CCF anatomy assets for web and local 3D rendering.

This script reads the official Allen 100-um annotation volume, extracts the
whole-brain and unilateral VPL surfaces, and stores both triangle meshes plus
bit-packed voxel masks.  Coordinates are converted to the documented generic
bregma QC frame.  A subject-specific registration must replace this generic
mapping for formal inference.
"""
from __future__ import annotations

import argparse
import base64
import json
from pathlib import Path

import numpy as np
import nrrd
from skimage.measure import marching_cubes
import trimesh


ALLEN_ANNOTATION_URL = "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_100.nrrd"
ALLEN_TEMPLATE_URL = "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_100.nrrd"
VPL_ID = 718
RESOLUTION_MM = 0.1
BREGMA_AP_MM = 5.4
BREGMA_ML_MM = 5.7


def vertices_to_bregma(vertices_ap_dv_ml_mm: np.ndarray) -> np.ndarray:
    v = np.asarray(vertices_ap_dv_ml_mm, dtype=float)
    return np.column_stack([
        BREGMA_AP_MM - v[:, 0],
        v[:, 2] - BREGMA_ML_MM,
        -v[:, 1],
    ])


def compact_mesh(mask: np.ndarray, face_count: int) -> dict:
    vertices, faces, _, _ = marching_cubes(mask.astype(np.uint8), 0.5, spacing=(RESOLUTION_MM,) * 3)
    mesh = trimesh.Trimesh(vertices=vertices_to_bregma(vertices), faces=faces, process=True)
    if len(mesh.faces) > face_count:
        mesh = mesh.simplify_quadric_decimation(face_count=face_count)
    mesh.remove_unreferenced_vertices()
    return {
        "vertices": np.round(mesh.vertices, 3).tolist(),
        "faces": mesh.faces.astype(int).tolist(),
    }


def pack_mask(mask: np.ndarray) -> str:
    packed = np.packbits(mask.ravel(order="C"), bitorder="little")
    return base64.b64encode(packed.tobytes()).decode("ascii")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--brain-faces", type=int, default=4200)
    parser.add_argument("--vpl-faces", type=int, default=1400)
    args = parser.parse_args()

    annotation, header = nrrd.read(args.annotation)
    if annotation.ndim != 3:
        raise SystemExit("Annotation must be 3-D")
    brain = annotation > 0
    ml_indices = np.indices(annotation.shape, sparse=True)[2]
    left_vpl = (annotation == VPL_ID) & (ml_indices < annotation.shape[2] // 2)
    if not np.any(left_vpl):
        raise SystemExit("VPL structure id 718 not present in annotation")

    left_indices = np.argwhere(left_vpl)
    ccf_mm = left_indices * RESOLUTION_MM
    vpl_points = vertices_to_bregma(ccf_mm)
    requested = np.asarray([-1.94, -1.75, -3.10])
    distances = np.linalg.norm(vpl_points - requested, axis=1)
    nearest = vpl_points[int(np.argmin(distances))]

    payload = {
        "schema": "allen-ccf-web-asset-v1",
        "source": {
            "atlas": "Allen Mouse CCFv3 2017 annotation",
            "annotation_url": ALLEN_ANNOTATION_URL,
            "template_url": ALLEN_TEMPLATE_URL,
            "resolution_um": 100,
            "vpl_structure_id": VPL_ID,
        },
        "coordinate_frame": {
            "name": "generic_bregma_qc_mm",
            "axis_order": ["AP", "ML", "DV"],
            "convention": "AP anterior+, ML right+, DV dorsal+",
            "ccf_bregma_mm_AP_DV_ML": [BREGMA_AP_MM, 0.0, BREGMA_ML_MM],
            "formal_use_warning": "Generic QC mapping only; use subject/template registration for formal FEM.",
        },
        "shape_AP_DV_ML": list(annotation.shape),
        "brain_mask_bitpack_little_base64": pack_mask(brain),
        "left_vpl_mask_bitpack_little_base64": pack_mask(left_vpl),
        "brain_mesh": compact_mesh(brain, args.brain_faces),
        "left_vpl_mesh": compact_mesh(left_vpl, args.vpl_faces),
        "requested_target_bregma_mm": requested.tolist(),
        "nearest_left_vpl_voxel_bregma_mm": np.round(nearest, 3).tolist(),
        "generic_mapping_target_to_vpl_distance_mm": float(distances.min()),
        "nrrd_header_space_directions": np.asarray(header.get("space directions", [])).tolist(),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Brain mesh: {len(payload['brain_mesh']['vertices'])} vertices, {len(payload['brain_mesh']['faces'])} faces")
    print(f"Left VPL mesh: {len(payload['left_vpl_mesh']['vertices'])} vertices, {len(payload['left_vpl_mesh']['faces'])} faces")
    print(f"Generic mapping target-to-VPL distance: {distances.min():.3f} mm")


if __name__ == "__main__":
    main()

