#!/usr/bin/env python
"""Build a conformal four-tissue ellipsoidal mouse-head demonstration mesh.

This mesh validates the numerical pipeline and web/FEM handoff.  It is not an
anatomical replacement for micro-CT/MRI/EFMouse geometry.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="outputs/demo_mouse_head.msh")
    parser.add_argument("--mesh-size-mm", type=float, default=0.35)
    args = parser.parse_args()

    try:
        import gmsh
    except ModuleNotFoundError as exc:
        raise SystemExit("gmsh 未安装：python -m pip install gmsh") from exc

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("demo_mouse_head_bregma")
        occ = gmsh.model.occ
        center = (-1.25, 0.0, -3.05)

        def ellipsoid(rx: float, ry: float, rz: float) -> int:
            tag = occ.addSphere(*center, 1.0)
            occ.dilate([(3, tag)], *center, rx, ry, rz)
            return tag

        brain = ellipsoid(4.45, 4.35, 2.65)
        csf = ellipsoid(4.60, 4.50, 2.78)
        skull = ellipsoid(4.85, 4.75, 2.96)
        scalp = ellipsoid(5.10, 5.00, 3.16)
        occ.fragment([(3, scalp)], [(3, skull), (3, csf), (3, brain)])
        occ.synchronize()

        volumes = gmsh.model.getEntities(3)
        if len(volumes) != 4:
            raise RuntimeError(f"Expected four partitioned tissue volumes, found {len(volumes)}")

        def outer_ap_span(entity: tuple[int, int]) -> float:
            box = gmsh.model.getBoundingBox(*entity)
            return box[3] - box[0]

        inner_to_outer = sorted(volumes, key=outer_ap_span)
        assignments = [
            (inner_to_outer[0][1], 4, "brain"),
            (inner_to_outer[1][1], 3, "csf"),
            (inner_to_outer[2][1], 2, "skull"),
            (inner_to_outer[3][1], 1, "scalp"),
        ]
        for entity_tag, physical_tag, name in assignments:
            gmsh.model.addPhysicalGroup(3, [entity_tag], physical_tag)
            gmsh.model.setPhysicalName(3, physical_tag, name)

        h = float(args.mesh_size_mm)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.72 * h)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1.25 * h)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(output))
    finally:
        gmsh.finalize()
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

