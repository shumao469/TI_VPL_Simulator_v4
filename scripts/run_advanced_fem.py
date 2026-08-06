#!/usr/bin/env python3
"""Run registered two-channel TI FEM on the official EFMouse tetrahedral mesh."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import meshio
import nrrd
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ti_vpl.envelope import carrier_rms, directional_envelope, maximum_orientation_envelope, vector_angle_deg
from ti_vpl.fem import boundary_faces_from_tets, read_gmsh_tetra, select_surface_patch, solve_current_pair, triangle_geometry
from ti_vpl.io_utils import element_centroids, roi_statistics, spherical_roi, write_json
from ti_vpl.registration import apply_affine, read_affine, validate_registration


TISSUE_IDS = {"gray": 1, "csf": 2, "bone": 3, "skin": 4, "eye": 5}


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def nearest_patch(mesh, boundary: np.ndarray, requested_fem_mm: np.ndarray, radius_mm: float) -> tuple[np.ndarray, float, np.ndarray]:
    centroids, _, _ = triangle_geometry(mesh.points_m, boundary)
    requested_m = np.asarray(requested_fem_mm) * 1e-3
    distances = np.linalg.norm(centroids - requested_m, axis=1)
    center_m = centroids[int(np.argmin(distances))]
    faces = select_surface_patch(mesh.points_m, boundary, center_m, radius_mm * 1e-3)
    return faces, float(distances.min() * 1e3), center_m * 1e3


def atlas_vpl_mask(centroids_fem_mm: np.ndarray, fem_from_bregma: np.ndarray, annotation_path: Path, origin: list[float]) -> np.ndarray:
    annotation, _ = nrrd.read(str(annotation_path), index_order="C")
    bregma = apply_affine(centroids_fem_mm, np.linalg.inv(fem_from_bregma))
    ap0, dv0, ml0 = [float(v) for v in origin]
    # Bundled annotation is 100 um. Array order is AP, DV, ML.
    ijk = np.rint(np.column_stack([(ap0 - bregma[:, 0]) / 0.1, -bregma[:, 2] / 0.1, (bregma[:, 1] + ml0) / 0.1])).astype(int)
    inside = np.all((ijk >= 0) & (ijk < np.asarray(annotation.shape)), axis=1)
    mask = np.zeros(len(ijk), dtype=bool)
    valid = ijk[inside]
    mask[inside] = annotation[valid[:, 0], valid[:, 1], valid[:, 2]] == 718
    return mask


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_advanced.yaml")
    parser.add_argument("--demo-unregistered", action="store_true", help="Explicitly allow identity transform; output is non-formal")
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    # Project-relative paths remain stable even when the web app writes a run-specific config.
    base = ROOT
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    mesh_path = resolve(base, cfg["anatomy"]["converted_mesh_vtu"])
    if not mesh_path.exists():
        raise SystemExit(f"Converted EFMouse mesh missing: {mesh_path}\nRun download_efmouse_assets.py then convert_efmouse_mesh.py")
    affine_path = resolve(base, cfg["coordinates"]["bregma_to_fem_affine_file"])
    qc_path = resolve(base, cfg["coordinates"]["registration_qc_file"])
    if not affine_path.exists() or not qc_path.exists():
        if not args.demo_unregistered:
            raise SystemExit("Formal FEM stopped: registration affine/QC missing. Run scripts/fit_registration.py first.")
        fem_from_bregma = np.eye(4)
        registration = {"passed": False, "mode": "DEMO_UNREGISTERED_IDENTITY", "rms_mm": None}
    else:
        fem_from_bregma = read_affine(affine_path)
        registration = json.loads(qc_path.read_text(encoding="utf-8"))
        validate_registration(registration, float(cfg["coordinates"]["maximum_registration_rms_mm"]))
        registration["mode"] = "REGISTERED"

    mesh = read_gmsh_tetra(mesh_path, cfg["anatomy"].get("coordinate_unit", "mm"))
    boundary = boundary_faces_from_tets(mesh.tets)
    stim = cfg["stimulation"]
    requested = []
    labels = []
    for channel in ("channel_1", "channel_2"):
        requested.extend([stim[channel]["source_bregma_mm"], stim[channel]["return_bregma_mm"]])
        labels.extend([f"{channel}_source", f"{channel}_return"])
    requested_fem = apply_affine(np.asarray(requested, dtype=float), fem_from_bregma)
    patches, projections, actual_centers = {}, {}, {}
    for label, center in zip(labels, requested_fem):
        faces, distance, actual = nearest_patch(mesh, boundary, center, float(stim["electrode_patch_radius_mm"]))
        patches[label], projections[label], actual_centers[label] = faces, distance, actual.tolist()
    max_projection = float(stim["maximum_surface_projection_mm"])
    if max(projections.values()) > max_projection and not args.demo_unregistered:
        raise SystemExit(f"Electrode-to-surface projection QC failed: {projections}; allowed <= {max_projection} mm")

    conductivities = {TISSUE_IDS[name]: float(value) for name, value in cfg["materials"]["conductivity_S_per_m"].items()}
    solver = cfg["solver"]
    common = dict(
        mesh=mesh,
        conductivity_by_tag=conductivities,
        gauge_node=int(solver["gauge_node"]),
        direct_max_nodes=int(solver["direct_maximum_nodes"]),
        rtol=float(solver["relative_tolerance"]),
        maxiter=int(solver["maximum_iterations"]),
    )
    current_1_A = float(stim.get("channel_1_current_mA", stim["current_per_channel_mA"])) * 1e-3
    current_2_A = float(stim.get("channel_2_current_mA", stim["current_per_channel_mA"])) * 1e-3
    ch1 = solve_current_pair(source_faces=patches["channel_1_source"], sink_faces=patches["channel_1_return"], current_A=current_1_A, **common)
    ch2 = solve_current_pair(source_faces=patches["channel_2_source"], sink_faces=patches["channel_2_return"], current_A=current_2_A, **common)
    allowed_imbalances = np.abs([current_1_A, current_2_A]) * float(solver["require_current_balance_relative"])
    imbalances = [float(ch1["boundary_current_imbalance_A"]), float(ch2["boundary_current_imbalance_A"])]
    if np.any(np.abs(imbalances) > allowed_imbalances):
        raise RuntimeError(f"Current balance QC failed: {imbalances} A")

    e1, e2 = np.asarray(ch1["E_V_per_m"]), np.asarray(ch2["E_V_per_m"])
    ti = maximum_orientation_envelope(e1, e2)
    centroids_m = element_centroids(mesh.points_m, mesh.tets)
    centroids_mm = centroids_m * 1e3
    target_cfg = cfg["coordinates"]["requested_target_mm"]
    target_bregma = np.asarray([target_cfg["AP"], target_cfg["ML"], target_cfg["DV"]], dtype=float)
    target_fem = apply_affine(target_bregma, fem_from_bregma)
    target_sphere = spherical_roi(centroids_m, target_fem, float(stim["target_sphere_radius_mm"])) & (mesh.cell_tags == 1)
    annotation_path = base / "assets/annotation_100.nrrd"
    vpl = atlas_vpl_mask(centroids_mm, fem_from_bregma, annotation_path, cfg["coordinates"]["generic_ccf_bregma_mm_AP_DV_ML"]) & (mesh.cell_tags == 1)
    gray = mesh.cell_tags == 1
    if not np.any(vpl):
        raise RuntimeError("Registered Allen VPL ROI contains no gray-matter tetrahedra; inspect registration and axis conventions")

    fields = {
        "tissue_id": [mesh.cell_tags.astype(np.int16)],
        "gray_matter": [gray.astype(np.uint8)],
        "Allen_VPL": [vpl.astype(np.uint8)],
        "target_sphere": [target_sphere.astype(np.uint8)],
        "E1_vector_V_per_m": [e1], "E2_vector_V_per_m": [e2],
        "E1_magnitude_V_per_m": [np.linalg.norm(e1, axis=1)], "E2_magnitude_V_per_m": [np.linalg.norm(e2, axis=1)],
        "TI_envelope_max_V_per_m": [ti],
        "TI_AP_V_per_m": [directional_envelope(e1, e2, fem_from_bregma[:3, 0])],
        "TI_ML_V_per_m": [directional_envelope(e1, e2, fem_from_bregma[:3, 1])],
        "TI_DV_V_per_m": [directional_envelope(e1, e2, fem_from_bregma[:3, 2])],
        "carrier_RMS_V_per_m": [carrier_rms(e1, e2)],
        "E1_E2_angle_deg": [vector_angle_deg(e1, e2)],
    }
    outdir = resolve(base, cfg["outputs"]["directory"])
    outdir.mkdir(parents=True, exist_ok=True)
    result_path = outdir / "ti_vpl_efmouse_registered.vtu"
    result = meshio.Mesh(points=mesh.points_m * 1e3, cells=[("tetra", mesh.tets)],
                         point_data={"phi1_V": ch1["phi_V"], "phi2_V": ch2["phi_V"]}, cell_data=fields)
    meshio.write(result_path, result, binary=True)
    metrics = {
        "scientific_status": "registered_template_FEM" if registration["mode"] == "REGISTERED" else "DEMO_UNREGISTERED_NOT_FOR_INFERENCE",
        "registration": registration,
        "coordinate_convention": "[AP anterior+, ML right+, DV dorsal+] mm from bregma",
        "requested_target_bregma_mm": target_bregma.tolist(), "target_fem_mm": target_fem.tolist(),
        "frequencies_Hz": [float(stim["f1_Hz"]), float(stim["f2_Hz"])],
        "beat_frequency_Hz": abs(float(stim["f2_Hz"]) - float(stim["f1_Hz"])),
        "channel_currents_mA": [current_1_A * 1e3, current_2_A * 1e3],
        "frequency_interpretation": "In this resistive electroquasistatic model, delta-f controls temporal beat rate; nonzero delta-f does not rescale the spatial envelope.",
        "electrode_surface_projection_mm": projections, "actual_surface_centers_fem_mm": actual_centers,
        "current_balance_A": imbalances,
        "linear_relative_residual": [float(ch1["linear_relative_residual"]), float(ch2["linear_relative_residual"])],
        "Allen_VPL": roi_statistics(ti, ch1["volumes_m3"], vpl),
        "target_sphere": roi_statistics(ti, ch1["volumes_m3"], target_sphere),
        "gray_matter": roi_statistics(ti, ch1["volumes_m3"], gray),
        "warning": "Registered template model, not an individualized animal dose or efficacy prediction.",
    }
    metrics["focality_VPL_median_over_gray_P95"] = metrics["Allen_VPL"]["median"] / max(metrics["gray_matter"]["p95"], 1e-30)
    write_json(outdir / "ti_vpl_efmouse_registered_metrics.json", metrics)
    print(f"Wrote {result_path}")
    print(f"Wrote {outdir / 'ti_vpl_efmouse_registered_metrics.json'}")


if __name__ == "__main__":
    main()
