#!/usr/bin/env python3
"""Convert EFMouse MATLAB 7.3/HDF5 node+elem arrays to a tagged VTU mesh."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import meshio
import numpy as np


def find_dataset(handle: h5py.File, basename: str) -> h5py.Dataset:
    found: list[h5py.Dataset] = []

    def visitor(name: str, obj) -> None:
        if isinstance(obj, h5py.Dataset) and name.rsplit("/", 1)[-1].lower() == basename.lower():
            found.append(obj)

    handle.visititems(visitor)
    if len(found) != 1:
        paths = [item.name for item in found]
        raise KeyError(f"Expected one dataset named {basename!r}, found {paths}")
    return found[0]


def normalize_node_elem(node: np.ndarray, elem: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    node = np.asarray(node)
    elem = np.asarray(elem)
    while node.ndim > 2 and 1 in node.shape:
        node = np.squeeze(node, axis=node.shape.index(1))
    while elem.ndim > 2 and 1 in elem.shape:
        elem = np.squeeze(elem, axis=elem.shape.index(1))
    if node.ndim != 2 or elem.ndim != 2:
        raise ValueError(f"node and elem must be 2D; got {node.shape}, {elem.shape}")
    if node.shape[1] == 3:
        points = node
    elif node.shape[0] == 3:
        points = node.T
    else:
        raise ValueError(f"Cannot identify node coordinate axis in shape {node.shape}")
    if elem.shape[1] == 5:
        records = elem
    elif elem.shape[0] == 5:
        records = elem.T
    else:
        raise ValueError(f"Cannot identify elem record axis in shape {elem.shape}")
    tets = np.rint(records[:, :4]).astype(np.int64)
    labels = np.rint(records[:, 4]).astype(np.int32)
    if tets.min() == 1:
        tets -= 1
    if tets.min() < 0 or tets.max() >= len(points):
        raise ValueError(f"Connectivity range [{tets.min()}, {tets.max()}] invalid for {len(points)} points")
    if not np.all(np.isfinite(points)):
        raise ValueError("Mesh contains non-finite coordinates")
    if not set(np.unique(labels)).issubset({1, 2, 3, 4, 5}):
        raise ValueError(f"Unexpected EFMouse tissue IDs: {np.unique(labels).tolist()}")
    p = points[tets]
    signed6 = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
    bad = np.count_nonzero(np.abs(signed6) <= 1e-15)
    if bad:
        raise ValueError(f"Mesh contains {bad} degenerate tetrahedra")
    return points.astype(np.float64), tets, labels


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/efmouse/Alekseichuk_mesh_clean.mat")
    parser.add_argument("--output", default="data/efmouse/Alekseichuk_mesh_clean.vtu")
    parser.add_argument("--metadata", default="data/efmouse/Alekseichuk_mesh_clean_qc.json")
    args = parser.parse_args()
    with h5py.File(args.input, "r") as handle:
        points, tets, labels = normalize_node_elem(find_dataset(handle, "node")[()], find_dataset(handle, "elem")[()])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    meshio.write_points_cells(output, points, [("tetra", tets)], cell_data={"tissue_id": [labels]}, binary=True)
    payload = {
        "source": str(Path(args.input).resolve()),
        "output": str(output.resolve()),
        "coordinate_unit": "mm (as distributed by EFMouse)",
        "points": len(points),
        "tetrahedra": len(tets),
        "bounds_mm": [points.min(axis=0).tolist(), points.max(axis=0).tolist()],
        "tissue_counts": {str(int(tag)): int(np.count_nonzero(labels == tag)) for tag in np.unique(labels)},
    }
    Path(args.metadata).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
