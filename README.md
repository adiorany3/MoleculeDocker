# Molecul Docking (Streamlit Demo)

This repository contains a demo Streamlit application for molecular docking. It supports two modes:

- AutoDock Vina backend (if `vina` binary is available and you provide receptor/ligand in PDBQT format)
- Mock docking backend using RDKit (generates conformers and scores them using UFF energies)

This is intended as a prototype and educational tool, not a production-ready docking server.

## Quickstart (macOS / Linux)

Install dependencies. Using conda for RDKit/Open Babel is recommended (most reliable):

```bash
# Recommended: create a conda env
conda create -n molecul python=3.10 -y
conda activate molecul
conda install -c conda-forge rdkit openbabel autodock-vina -y
pip install -r requirements.txt

Alternative: install all conda requirements from `requirements-conda.txt`:

```bash
conda install -c conda-forge --file requirements-conda.txt
```
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

## Results & Insight

- The app now shows a **Summary** panel with a concise analysis of the docking run: count, best score (kcal/mol), average score, score range and standard deviation. The best (lowest) score is marked with a trophy icon 🏆.
- The summary includes a simple heuristic interpretation (🔬/⭐/⚠️/❌) to quickly flag promising vs weak binders. Lower (more negative) scores indicate stronger predicted binding.
- Download buttons include format labels and a download icon (⬇️) for quick recognition.
 - A sidebar control allows choosing the image width for pose thumbnails.
 - Each pose shows a score-category badge (🔥/⭐/⚠️/❌) based on the score thresholds and the top pose is marked with a trophy 🏆.
 - Summary metrics (Best, Average, Std Dev) are shown as visual cards for quick inspection.
 - The app shows a **Best pose interpretation** panel for the top pose, including RDKit-based molecular descriptors (MW, H-bond donors/acceptors, TPSA, rotatable bonds, ligand efficiency) when RDKit is available — otherwise it shows PDB-based atom counts and a size estimate.

How to improve docking results — practical tips:

- **Prepare receptor & ligand carefully**: Remove crystallographic waters from the receptor, add hydrogens, assign protonation states (pH 7.4) and generate PDBQT files — these preprocessing steps reduce artifacts.
- **Grid box selection**: Use `Compute grid from receptor` or set the center and size to cover the known binding pocket. Too small boxes truncate poses; too large boxes increase search space and may reduce effective sampling.
- **Increase sampling**: Raise `exhaustiveness` (Vina) and `Number of poses` to explore more conformational space. Typical `exhaustiveness` values: 8–32. Double-check run times.
- **Use multiple scoring functions / consensus**: Re-score top poses with alternative scoring functions or tools to reduce false positives and increase confidence.
- **Clustering + consensus**: If pose variability is large, perform pose clustering (RMSD-based) and use consensus pose or representative centroid for analysis.
- **Apply medicinal chemistry filters**: Check ligand properties (MW, TPSA, HBD/HBA, rotatable bonds, ligand efficiency) before investing in experimental validation.
- **Check contacts**: If the app reports few contacts between the ligand and receptor, adjust grid or ligand placement — adding anchor points or constraints can help direct docking.

Example: for a small ligand such as ethanol (SMILES: `CCO`) a strong score like `-8.0` yields ligand efficiency (LE) ≈ -2.67 (kcal/mol per heavy atom), which falls in the "Excellent" range and is labeled with 🔥.

Quick snippet (Python) to compute a sample interpretation in a local REPL:
```python
from rdkit import Chem
from docking.insights import compute_pose_interpretation
mol = Chem.MolFromSmiles('CCO')
mol = Chem.AddHs(mol)
pose = {'mol': mol, 'score': -8.0}
interp = compute_pose_interpretation(pose)
print(interp['text'])
```

Extra features:
- **Convert to PDBQT**: use the Open Babel `obabel` CLI to convert uploaded PDB/MOL/SDF files to PDBQT (sidebar Convert to PDBQT button). If `obabel` isn't found, the app will show an error and fallback to mock docking.
- **3D Viewer**: enable `Show 3D viewer (py3Dmol)` in the sidebar to visualize poses in the browser when py3Dmol is installed. Without RDKit the mock PDB pose still displays as a 3D model.
 - **Compute grid**: use `Compute grid from receptor` in the sidebar to compute the docking box center and size from a receptor PDB. This automatically fills the `Center X/Y/Z` and `Size X/Y/Z` fields. (The calculation adds a small margin to the bounding box.)

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

Compiled by Galuh Adi Insani
