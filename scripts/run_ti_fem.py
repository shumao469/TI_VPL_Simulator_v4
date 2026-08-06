#!/usr/bin/env python
"""Solve two carrier channels and export TI-VPL FEM fields."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import meshio
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ti_vpl.coordinates import BregmaPoint, crossed_montage
from ti_vpl.envelope import carrier_rms, directional_envelope, maximum_orientation_envelope, vector_angle_deg
from ti_vpl.fem import boundary_faces_from_tets, read_gmsh_tetra, select_surface_patch, solve_current_pair, triangle_geometry
from ti_vpl.io_utils import element_centroids, roi_statistics, spherical_roi, write_json


def nearest_surface_center(points_m: np.ndarray, faces: np.ndarray, requested_mm: np.ndarray) -> tuple[np.ndarray, float]:
    centroids, _, _ = triangle_geometry(points_m, faces)
    requested_m = np.asarray(requested_mm, dtype=float) * 1e-3
    distances = np.linalg.norm(centroids - requested_m[None, :], axis=1)
    index = int(np.argmin(distances))
    return centroids[index], float(distances[index] * 1e3)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_vpl.yaml")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    base = config_path.parent
    mesh_path = Path(cfg["mesh"]["path"])
    if not mesh_path.is_absolute():
        mesh_path = base / mesh_path
    if not mesh_path.exists():
        raise SystemExit(
            f"Mesh not found: {mesh_path}\n"
            f"Run: python scripts/build_demo_mesh.py --output {mesh_path}"
        )

    t = cfg["coordinates"]["target_vpl_mm"]
    target = BregmaPoint(float(t["AP"]), float(t["ML"]), float(t["DV"]))
    stim = cfg["stimulation"]
    montage_cfg = stim["montage"]
    montage = crossed_montage(
        target,
        float(montage_cfg["half_spacing_mm"]),
        float(montage_cfg["electrode_dv_mm"]),
        float(montage_cfg.get("rotation_deg", 0.0)),
    )

    mesh = read_gmsh_tetra(mesh_path, cfg["mesh"].get("coordinate_unit", "mm"))
    boundary = boundary_faces_from_tets(mesh.tets)
    projected: dict[str, list[np.ndarray]] = {}
    projection_distances: dict[str, list[float]] = {}
    patch_faces: dict[str, list[np.ndarray]] = {}
    radius_m = float(montage_cfg["patch_radius_mm"]) * 1e-3
    for channel_name, pair in montage.items():
        projected[channel_name] = []
        projection_distances[channel_name] = []
        patch_faces[channel_name] = []
        for point in pair:
            centre, distance_mm = nearest_surface_center(mesh.points_m, boundary, point.as_array())
            projected[channel_name].append(centre)
            projection_distances[channel_name].append(distance_mm)
            patch_faces[channel_name].append(select_surface_patch(mesh.points_m, boundary, centre, radius_m))

    names = cfg["materials"]["physical_tag_to_name"]
    conductivity_named = cfg["materials"]["conductivity_S_per_m"]
    conductivity = {int(tag): float(conductivity_named[name]) for tag, name in names.items()}
    solver = cfg["solver"]
    solve_kwargs = dict(
        conductivity_by_tag=conductivity,
        current_A=float(stim["normalization_current_mA"]) * 1e-3,
        gauge_node=int(solver["gauge_node"]),
        direct_max_nodes=int(solver["direct_max_nodes"]),
        rtol=float(solver["iterative_rtol"]),
        maxiter=int(solver["iterative_maxiter"]),
    )
    ch1 = solve_current_pair(
        mesh,
        source_faces=patch_faces["channel_1"][0],
        sink_faces=patch_faces["channel_1"][1],
        **solve_kwargs,
    )
    ch2 = solve_current_pair(
        mesh,
        source_faces=patch_faces["channel_2"][0],
        sink_faces=patch_faces["channel_2"][1],
        **solve_kwargs,
    )

    e1 = np.asarray(ch1["E_V_per_m"], dtype=float)
    e2 = np.asarray(ch2["E_V_per_m"], dtype=float)
    mti = maximum_orientation_envelope(e1, e2)
    centroids = element_centroids(mesh.points_m, mesh.tets)
    roi = spherical_roi(centroids, target.as_array(), float(cfg["coordinates"]["target_radius_mm"]))
    brain = mesh.cell_tags == 4
    roi &= brain
    volumes = ch1["volumes_m3"]

    outdir = Path(cfg["outputs"]["directory"])
    if not outdir.is_absolute():
        outdir = base / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    cell_data = {
        "material_tag": [mesh.cell_tags.astype(np.int32)],
        "brain": [brain.astype(np.uint8)],
        "VPL_sphere": [roi.astype(np.uint8)],
        "E1_magnitude_V_per_m": [np.linalg.norm(e1, axis=1)],
        "E2_magnitude_V_per_m": [np.linalg.norm(e2, axis=1)],
        "TI_envelope_max_V_per_m": [mti],
        "TI_AP_V_per_m": [directional_envelope(e1, e2, np.asarray([1.0, 0.0, 0.0]))],
        "TI_ML_V_per_m": [directional_envelope(e1, e2, np.asarray([0.0, 1.0, 0.0]))],
        "TI_DV_V_per_m": [directional_envelope(e1, e2, np.asarray([0.0, 0.0, 1.0]))],
        "carrier_RMS_V_per_m": [carrier_rms(e1, e2)],
        "E1_E2_angle_deg": [vector_angle_deg(e1, e2)],
        "E1_vector_V_per_m": [e1],
        "E2_vector_V_per_m": [e2],
    }
    result_mesh = meshio.Mesh(
        points=mesh.points_m * 1e3,
        cells=[("tetra", mesh.tets)],
        point_data={"phi1_V": ch1["phi_V"], "phi2_V": ch2["phi_V"]},
        cell_data=cell_data,
    )
    result_path = outdir / "ti_vpl_fields.vtu"
    meshio.write(result_path, result_mesh)

    brain_stats = roi_statistics(mti, volumes, brain)
    vpl_stats = roi_statistics(mti, volumes, roi)
    metrics = {
        "coordinate_frame": "bregma_relative_mm; AP anterior+, ML right+, DV dorsal+",
        "target_vpl_mm": target.as_dict(),
        "target_radius_mm": float(cfg["coordinates"]["target_radius_mm"]),
        "frequencies_Hz": [float(stim["f1_Hz"]), float(stim["f2_Hz"])],
        "beat_frequency_Hz": abs(float(stim["f2_Hz"]) - float(stim["f1_Hz"])),
        "normalization_current_mA": float(stim["normalization_current_mA"]),
        "vpl": vpl_stats,
        "brain": brain_stats,
        "focality_vpl_median_over_brain_p95": float(vpl_stats["median"] / max(float(brain_stats["p95"]), 1e-30)),
        "electrode_surface_projection_distance_mm": projection_distances,
        "warning": "Demonstration ellipsoid FEM; replace geometry, contacts and registered ROI for publication or experimental planning.",
    }
    write_json(outdir / "ti_vpl_metrics.json", metrics)
    print(f"Wrote {result_path}")
    print(f"Wrote {outdir / 'ti_vpl_metrics.json'}")
    print(f"VPL median mTI = {float(vpl_stats['median']):.6g} V/m")


if __name__ == "__main__":
    main()
