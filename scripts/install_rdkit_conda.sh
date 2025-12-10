#!/usr/bin/env bash
set -euo pipefail

# Simple installer for the Molecul demo using conda (miniforge recommended)
# Usage: bash scripts/install_rdkit_conda.sh

ENV_NAME=${1:-molecul}
PY_VERSION=${2:-3.10}

if ! command -v conda >/dev/null 2>&1; then
  echo "Conda not found. Please install Miniforge (https://github.com/conda-forge/miniforge) or Anaconda first."
  exit 1
fi

echo "Creating conda env: $ENV_NAME (python $PY_VERSION)"
conda create -n "$ENV_NAME" python="$PY_VERSION" -y

echo "Activating env $ENV_NAME"
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "Installing RDKit/Open Babel/Vina from conda-forge and other deps"
conda install -c conda-forge --file requirements-conda.txt -y
pip install -r requirements.txt

echo "Done. Activate with: conda activate $ENV_NAME"
