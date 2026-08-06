# EFMouse data directory

Run `python scripts/download_efmouse_assets.py` and then
`python scripts/convert_efmouse_mesh.py` to populate this directory.

The source MAT and converted VTU are intentionally excluded from the release
ZIP because the official source mesh is about 94 MB. The included manifest and
QC JSON records the exact official file validated while preparing Xu Lab v4.
