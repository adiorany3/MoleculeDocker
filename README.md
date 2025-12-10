# Molecul Docking (Streamlit Demo)

This repository contains a demo Streamlit application for molecular docking. It supports two modes:

- AutoDock Vina backend (if `vina` binary is available and you provide receptor/ligand in PDBQT format)
- Mock docking backend using RDKit (generates conformers and scores them using UFF energies)

This is intended as a prototype and educational tool, not a production-ready docking server.

## Quickstart (macOS / Linux)

Install dependencies. Using conda for RDKit is recommended (most reliable):

```bash
# Recommended: create a conda env
conda create -n molecul python=3.10 -y
conda activate molecul
conda install -c conda-forge rdkit -y
pip install -r requirements.txt
```

Or with pip only (RDKit via rdkit-pypi may not include all features or be slow to install):

```bash
pip install -r requirements.txt

Troubleshooting note: On macOS (especially Apple Silicon / M1/M2) and some Linux platforms, `rdkit-pypi` is not published for all Python versions/architectures, hence you may encounter:

```
ERROR: Could not find a version that satisfies the requirement rdkit-pypi
ERROR: No matching distribution found for rdkit-pypi
```

If you see that, use conda (miniforge/conda-forge) instead — it provides pre-built RDKit packages for many platforms:

```bash
# miniforge recommended for Apple Silicon
conda create -n molecul python=3.10 -y
conda activate molecul
conda install -c conda-forge rdkit -y
pip install -r requirements.txt
```

If you prefer containers, use the provided `Dockerfile` which installs RDKit via conda-forge and sets up the environment.
```

Note: If you cannot install RDKit, the app includes a simple RDKit-free mock docking fallback that generates placeholder PDB poses and scores so the UI remains functional.

Install AutoDock Vina (optional): download and place `vina` on PATH.
Install Open Babel (optional; required for Convert to PDBQT):

```bash
# Conda (recommended):
conda activate molecul
conda install -c conda-forge openbabel -y

# macOS (Homebrew):
brew install open-babel

# Debian/Ubuntu:
sudo apt-get install openbabel -y
```

Verify RDKit installation (optional):

```bash
python -c "from rdkit import Chem; print(Chem.__name__)"
```

Run the app:

```bash
streamlit run app.py
```

Run tests (local):

```bash
bash scripts/run_tests.sh
```

CI: a GitHub Actions workflow is included at `.github/workflows/ci.yml` that sets up RDKit via conda and runs `pytest`.

## Usage

- Upload receptor and ligand files on the left-hand sidebar.
- If Vina is available and both files are PDBQT, Vina will be used.
- Otherwise, the app uses a mock docking algorithm with RDKit.
- It displays poses and allows downloading SDF files for poses.

Extra features:
- **Convert to PDBQT**: use the Open Babel `obabel` CLI to convert uploaded PDB/MOL/SDF files to PDBQT (sidebar Convert to PDBQT button). If `obabel` isn't found, the app will show an error and fallback to mock docking.
- **3D Viewer**: enable `Show 3D viewer (py3Dmol)` in the sidebar to visualize poses in the browser when py3Dmol is installed. Without RDKit the mock PDB pose still displays as a 3D model.

## Notes and Limitations

- If you plan to do real docking, prepare receptor/ligand files in PDBQT using AutoDock Tools (MGLTools) or Open Babel.
- This demo uses RDKit to create conformers and compute UFF energies — these are not directly comparable to Vina scores.
- To use Vina in this app, ensure you upload `PDBQT` files.

## Files

- `app.py` — Streamlit app
- `docking/vina_runner.py` — wrapper for calling `vina`
- `docking/mock_dock.py` — RDKit-based mock docking
- `docking/utils.py` — RDKit utilities

## Next Steps

- Add 3D visualization using py3Dmol or nglview via Streamlit components.
- Add preprocessing and conversions to/from PDBQT using Open Babel or MGLTools.
- Add more advanced scoring and pose clustering.

If you want, I can now:
- Add 3D viewers
- Add a converter from PDB/PDBQT using Open Babel if available
- Add Dockerfile and CI

Tell me which next step you'd like me to take.
