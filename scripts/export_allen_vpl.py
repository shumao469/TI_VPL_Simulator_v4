#!/usr/bin/env python3
"""Download official Allen CCF NRRDs and export VPL without AllenSDK."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import urllib.request

import numpy as np
import nrrd

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION_URL = "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_{resolution}.nrrd"
TEMPLATE_URL = "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_{resolution}.nrrd"
GRAPH_URL = "https://api.brain-map.org/api/v2/structure_graph_download/1.json"


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 10_000:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "ti-vpl-fem/3.0"})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.replace(path)


def find_structure(graph: dict, acronym: str) -> dict:
    stack = list(graph["msg"])
    while stack:
        item = stack.pop()
        if item.get("acronym", "").upper() == acronym.upper():
            return item
        stack.extend(item.get("children", []))
    raise KeyError(f"Allen structure acronym not found: {acronym}")


def descendant_ids(item: dict) -> list[int]:
    return [int(item["id"]), *[child_id for child in item.get("children", []) for child_id in descendant_ids(child)]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="allen_ccf")
    parser.add_argument("--resolution-um", type=int, choices=[10, 25, 50, 100], default=100)
    parser.add_argument("--hemisphere", choices=["both", "left", "right"], default="left")
    parser.add_argument("--acronym", default="VPL")
    parser.add_argument("--target-ap-mm", type=float, default=-1.94)
    parser.add_argument("--target-ml-mm", type=float, default=-1.75)
    parser.add_argument("--target-dv-mm", type=float, default=-3.10)
    args = parser.parse_args()
    out = Path(args.output_dir)
    source = out / "source"
    source.mkdir(parents=True, exist_ok=True)

    bundled_annotation = ROOT / "assets/annotation_100.nrrd"
    bundled_template = ROOT / "assets/average_template_100.nrrd"
    if args.resolution_um == 100 and bundled_annotation.exists() and bundled_template.exists():
        annotation_path, source_template = bundled_annotation, bundled_template
    else:
        annotation_path = source / f"annotation_{args.resolution_um}.nrrd"
        source_template = source / f"average_template_{args.resolution_um}.nrrd"
        download(ANNOTATION_URL.format(resolution=args.resolution_um), annotation_path)
        download(TEMPLATE_URL.format(resolution=args.resolution_um), source_template)
    graph_path = source / "structure_graph_1.json"
    download(GRAPH_URL, graph_path)
    graph = json.loads(graph_path.read_text(encoding="utf-8"))
    structure = find_structure(graph, args.acronym)
    ids = descendant_ids(structure)

    annotation, header = nrrd.read(str(annotation_path), index_order="C")
    mask = np.isin(annotation, ids)
    mid = mask.shape[2] // 2  # CCF array order AP, DV, ML; first ML half is left.
    if args.hemisphere == "left":
        mask[:, :, mid:] = False
    elif args.hemisphere == "right":
        mask[:, :, :mid] = False
    mask = mask.astype(np.uint8)

    mask_path = out / f"{args.acronym}_{args.hemisphere}_{args.resolution_um}um.nrrd"
    template_path = out / f"template_{args.resolution_um}um.nrrd"
    nrrd.write(str(mask_path), mask, header=header, index_order="C")
    shutil.copy2(source_template, template_path)
    voxel_mm = args.resolution_um / 1000.0
    coordinates = np.argwhere(mask)
    centroid = coordinates.mean(axis=0) * voxel_mm if len(coordinates) else np.full(3, np.nan)
    stl_path = out / f"{args.acronym}_{args.hemisphere}_{args.resolution_um}um_visualization.stl"
    if len(coordinates):
        from skimage.measure import marching_cubes
        import trimesh

        vertices, faces, _, _ = marching_cubes(mask, 0.5, spacing=(voxel_mm,) * 3)
        trimesh.Trimesh(vertices=vertices, faces=faces, process=False).export(stl_path)
    info = {
        "source": "Allen CCFv3 2017 official current-release NRRD",
        "source_urls": [ANNOTATION_URL.format(resolution=args.resolution_um), TEMPLATE_URL.format(resolution=args.resolution_um), GRAPH_URL],
        "acronym": args.acronym, "structure_id": int(structure["id"]), "included_structure_ids": ids,
        "hemisphere": args.hemisphere, "resolution_um": args.resolution_um,
        "ccf_array_axis_order": ["AP_anterior_to_posterior", "DV_dorsal_to_ventral", "ML_left_to_right"],
        "mask_shape": list(mask.shape), "mask_voxels": int(mask.sum()), "volume_mm3": float(mask.sum() * voxel_mm**3),
        "centroid_ccf_array_mm_AP_DV_ML": centroid.tolist(),
        "requested_stereotaxic_target_bregma_mm_AP_ML_DV": [args.target_ap_mm, args.target_ml_mm, args.target_dv_mm],
        "registration_warning": "Bregma coordinates are not CCF voxels; apply a validated registration before combining the ROI with FEM.",
        "dependency_note": "No AllenSDK; PyPI package pynrrd imports as nrrd.",
    }
    (out / "VPL_info.json").write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(info, ensure_ascii=False, indent=2))
    print(f"Wrote {mask_path}\nWrote {template_path}\nWrote {stl_path}")


if __name__ == "__main__":
    main()
