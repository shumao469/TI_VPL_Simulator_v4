#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="ti_vpl_fem311"
ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

if command -v mamba >/dev/null 2>&1; then
  CONDA_CMD="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_CMD="conda"
else
  echo "ERROR: mamba/conda not found. Install Miniforge, then rerun." >&2
  exit 2
fi

if "$CONDA_CMD" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Environment $ENV_NAME already exists; updating it in place."
  "$CONDA_CMD" env update -n "$ENV_NAME" -f "$ROOT_DIR/environment.yml" --prune
else
  "$CONDA_CMD" env create -f "$ROOT_DIR/environment.yml"
fi

"$CONDA_CMD" run -n "$ENV_NAME" python -m pip install --no-deps -e "$ROOT_DIR"

cd "$ROOT_DIR"
"$CONDA_CMD" run -n "$ENV_NAME" python - <<'PY'
import sys
import numpy, scipy, h5py, nrrd, meshio, pyvista, streamlit
from brainglobe_atlasapi import BrainGlobeAtlas
print("Environment OK")
print("Python:", sys.version)
print("NumPy:", numpy.__version__)
print("SciPy:", scipy.__version__)
print("PyVista:", pyvista.__version__)
print("Streamlit:", streamlit.__version__)
print("BrainGlobe Atlas API import: OK")
PY

echo
echo "Activate with: conda activate $ENV_NAME"
