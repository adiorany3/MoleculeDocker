#!/usr/bin/env bash
set -euo pipefail

# Simple installer for Open Babel using conda (miniforge recommended)
# Usage: bash scripts/install_openbabel_conda.sh

ENV_NAME=${1:-molecul}

if ! command -v conda >/dev/null 2>&1; then
  echo "Conda not found. Please install Miniforge (https://github.com/conda-forge/miniforge) or Anaconda first."
  exit 1
fi

# Activate conda
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME" || {
  echo "Environment $ENV_NAME doesn't exist. Create it first or omit ENV_NAME to use default 'molecul'"
  exit 1
}

echo "Installing Open Babel via conda-forge..."
conda install -c conda-forge openbabel -y

echo "Open Babel installed. Verify with: obabel -V or python -c 'from openbabel import pybel; print(pybel)'
"
