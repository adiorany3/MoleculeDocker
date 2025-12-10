#!/usr/bin/env bash
set -euo pipefail

# Script to run the project tests.
# It tries to run pytest; if pytest isn't available, it offers to install via pip or conda.

AUTO_INSTALL=${1:-}
export PYTHONPATH="$(pwd)"
if command -v pytest >/dev/null 2>&1; then
  echo "Running pytest..."
  pytest -q
  exit 0
fi

# No pytest: try python -m pytest
if python -m pytest -q; then
  exit 0
fi

# If we reach here, pytest is likely not installed.

echo "pytest is not available in your environment."
if [[ "$AUTO_INSTALL" == "--auto" ]]; then
  ans=Y
else
  read -p "Attempt to install pytest via pip now? [y/N] " ans
fi
if [[ "$ans" =~ ^[Yy]$ ]]; then
  echo "Installing pytest via pip..."
  python -m pip install --upgrade pip
  python -m pip install pytest
  pytest -q
  exit 0
fi

if [[ "$AUTO_INSTALL" == "--auto" ]]; then
  ans=Y
else
  read -p "Attempt to install pytest via conda (recommended if using conda) [y/N] " ans
fi
if [[ "$ans" =~ ^[Yy]$ ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "conda not available. Please install Miniforge or Anaconda, then run the script again.";
    exit 1
  fi
  echo "Installing pytest via conda..."
  conda install -c conda-forge pytest -y
  pytest -q
  exit 0
fi

echo "Tests not run. Install pytest and re-run: pip install pytest or conda install -c conda-forge pytest";
exit 1
