#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python -m streamlit run app_advanced.py --server.address 0.0.0.0 --server.port 8504
