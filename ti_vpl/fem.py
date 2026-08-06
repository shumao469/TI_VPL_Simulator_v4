"""Transparent P1 tetrahedral FEM solver for electro-quasistatic stimulation.

For each carrier channel the conductive model is
    div(sigma grad(phi)) = 0,  E = -grad(phi).
The source and return contacts use uniform Neumann current density.  This
minimal solver is intended for reproducible mouse-scale templates and method
development; large production meshes should use DOLFINx/PETSc or another
validated parallel solver.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import warnings

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import LinearOperator, cg, gmres, spsolve


@dataclass
class TetMesh:
    points_m: np.ndarray
    tets: np.ndarray
    cell_tags: np.ndarray


def read_gmsh_tetra(path: str | Path, coordinate_unit: str = "mm") -> TetMesh:
    import meshio

    mesh = meshio.read(str(path))
    physical = mesh.cell_data.get("tissue_id", mesh.cell_data.get("gmsh:physical"))
    tets, tags = [], []
    for block_i, block in enumerate(mesh.cells):
        if block.type not in {"tetra", "tetra10"}:
            continue
        data = np.asarray(block.data, dtype=np.int64)
        if block.type == "tetra10":
            warnings.warn("tetra10 found; using only the four corner nodes")
            data = data[:, :4]
        tets.append(data)
        tags.append(
            np.ones(len(data), dtype=np.int32)
            if physical is None
            else np.asarray(physical[block_i], dtype=np.int32)
        )
    if not tets:
        raise ValueError(f"No tetrahedra found in {path}")
    scale = {"m": 1.0, "mm": 1e-3, "um": 1e-6}[coordinate_unit]
    return TetMesh(
        points_m=np.asarray(mesh.points[:, :3], dtype=float) * scale,
        tets=np.vstack(tets),
        cell_tags=np.concatenate(tags),
    )


def tetra_geometry(points_m: np.ndarray, tets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = points_m[tets]
    d = np.empty((len(tets), 4, 4), dtype=float)
    d[:, :, 0] = 1.0
    d[:, :, 1:] = x
    inv_d = np.linalg.inv(d)
    grads = np.transpose(inv_d[:, 1:, :], (0, 2, 1))
    edges = x[:, 1:, :] - x[:, :1, :]
    volumes = np.abs(np.linalg.det(edges)) / 6.0
    if np.any(volumes <= 0):
        raise ValueError("Degenerate tetrahedra detected")
    return volumes, grads


def boundary_faces_from_tets(tets: np.ndarray) -> np.ndarray:
    local = np.asarray([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int64)
    faces = tets[:, local].reshape(-1, 3)
    canonical = np.sort(faces, axis=1)
    _, inverse, counts = np.unique(canonical, axis=0, return_inverse=True, return_counts=True)
    return faces[counts[inverse] == 1]


def triangle_geometry(points_m: np.ndarray, faces: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = points_m[faces]
    cross = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    normals = cross / np.maximum(2.0 * areas[:, None], 1e-30)
    return p.mean(axis=1), areas, normals


def select_surface_patch(
    points_m: np.ndarray,
    boundary_faces: np.ndarray,
    center_m: Iterable[float],
    radius_m: float,
) -> np.ndarray:
    center = np.asarray(center_m, dtype=float)
    centroids, areas, _ = triangle_geometry(points_m, boundary_faces)
    distance = np.linalg.norm(centroids - center[None, :], axis=1)
    keep = distance <= float(radius_m)
    if not np.any(keep):
        nearest = float(distance.min())
        raise ValueError(
            f"No surface triangles within radius; nearest face is {nearest * 1e3:.3f} mm away. "
            "Move the electrode to the actual mesh surface or enlarge patch_radius_mm."
        )
    if float(areas[keep].sum()) <= 0:
        raise ValueError("Selected electrode patch has zero area")
    return boundary_faces[keep]


def conductivity_per_element(cell_tags: np.ndarray, values: dict[int, float | complex]) -> np.ndarray:
    missing = sorted(set(np.unique(cell_tags)) - set(values))
    if missing:
        raise KeyError(f"Missing conductivity/admittivity for tags {missing}")
    out = np.empty(len(cell_tags), dtype=np.complex128)
    for tag, value in values.items():
        out[cell_tags == int(tag)] = value
    return out.real if np.all(np.abs(out.imag) < 1e-30) else out


def assemble_stiffness(
    points_m: np.ndarray,
    tets: np.ndarray,
    gamma_element: np.ndarray,
) -> tuple[csr_matrix, np.ndarray, np.ndarray]:
    volumes, grads = tetra_geometry(points_m, tets)
    gg = np.einsum("eik,ejk->eij", grads, grads)
    ke = gamma_element[:, None, None] * volumes[:, None, None] * gg
    rows = np.repeat(tets, 4, axis=1).ravel()
    cols = np.tile(tets, (1, 4)).ravel()
    matrix = coo_matrix((ke.reshape(-1), (rows, cols)), shape=(len(points_m), len(points_m))).tocsr()
    matrix.sum_duplicates()
    return matrix, volumes, grads


def add_uniform_current_load(
    rhs: np.ndarray,
    points_m: np.ndarray,
    faces: np.ndarray,
    current_A_entering: float,
) -> None:
    _, areas, _ = triangle_geometry(points_m, faces)
    q = float(current_A_entering) / float(areas.sum())
    per_vertex = q * areas / 3.0
    for local_vertex in range(3):
        np.add.at(rhs, faces[:, local_vertex], per_vertex)


def apply_gauge(matrix: csr_matrix, rhs: np.ndarray, node: int = 0) -> tuple[csr_matrix, np.ndarray]:
    a = matrix.tolil(copy=True)
    a[node, :] = 0
    a[:, node] = 0
    a[node, node] = 1
    b = rhs.copy()
    b[node] = 0
    return a.tocsr(), b


def solve_linear_system(
    matrix: csr_matrix,
    rhs: np.ndarray,
    direct_max_nodes: int = 300_000,
    rtol: float = 1e-10,
    maxiter: int = 8000,
) -> np.ndarray:
    if matrix.shape[0] <= direct_max_nodes:
        return spsolve(matrix, rhs)
    try:
        import pyamg

        ml = pyamg.smoothed_aggregation_solver(matrix)
        preconditioner = ml.aspreconditioner()
    except (ImportError, TypeError, ValueError):
        diag = matrix.diagonal()
        inv_diag = np.zeros_like(diag)
        nz = np.abs(diag) > 1e-30
        inv_diag[nz] = 1.0 / diag[nz]
        preconditioner = LinearOperator(matrix.shape, matvec=lambda x: inv_diag * x)
    if np.iscomplexobj(matrix.data):
        solution, info = gmres(matrix, rhs, M=preconditioner, rtol=rtol, atol=0.0, maxiter=maxiter)
    else:
        solution, info = cg(matrix, rhs, M=preconditioner, rtol=rtol, atol=0.0, maxiter=maxiter)
    if info != 0:
        raise RuntimeError(f"Iterative solver failed to converge; info={info}")
    return solution


def element_electric_field(phi: np.ndarray, tets: np.ndarray, grads: np.ndarray) -> np.ndarray:
    return -np.einsum("ei,eik->ek", phi[tets], grads)


def solve_current_pair(
    mesh: TetMesh,
    conductivity_by_tag: dict[int, float | complex],
    source_faces: np.ndarray,
    sink_faces: np.ndarray,
    current_A: float,
    gauge_node: int = 0,
    direct_max_nodes: int = 300_000,
    rtol: float = 1e-10,
    maxiter: int = 8000,
) -> dict[str, np.ndarray]:
    gamma = conductivity_per_element(mesh.cell_tags, conductivity_by_tag)
    matrix, volumes, grads = assemble_stiffness(mesh.points_m, mesh.tets, gamma)
    rhs = np.zeros(len(mesh.points_m), dtype=np.result_type(matrix.dtype, np.float64))
    add_uniform_current_load(rhs, mesh.points_m, source_faces, +current_A)
    add_uniform_current_load(rhs, mesh.points_m, sink_faces, -current_A)
    imbalance = float(abs(rhs.sum()))
    if imbalance > max(1e-14, abs(current_A) * 1e-8):
        warnings.warn(f"Discrete boundary current imbalance: {imbalance:.3e} A")
    matrix, rhs = apply_gauge(matrix, rhs, int(gauge_node))
    phi = solve_linear_system(matrix, rhs, direct_max_nodes, rtol, maxiter)
    return {
        "phi_V": phi,
        "E_V_per_m": element_electric_field(phi, mesh.tets, grads),
        "volumes_m3": volumes,
        "boundary_current_imbalance_A": np.asarray(imbalance),
        "linear_relative_residual": np.asarray(
            np.linalg.norm(matrix @ phi - rhs) / max(np.linalg.norm(rhs), 1e-30)
        ),
    }
