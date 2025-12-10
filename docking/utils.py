from io import BytesIO

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except Exception:
    Chem = None
    Draw = None
try:
    import py3Dmol
except Exception:
    py3Dmol = None

ALLOWED_EXT = [".pdb", ".pdbqt", ".mol", ".sdf"]


def allowed_file_ext(filename):
    import os
    return os.path.splitext(filename)[1].lower() in ALLOWED_EXT


def load_molecule_image(mol, confId=None, size=(300, 300)):
    if Draw is None:
        raise ImportError("RDKit not installed; cannot draw molecule")
    if confId is not None:
        img = Draw.MolToImage(mol, size=size, kekulize=True, confId=confId)
    else:
        img = Draw.MolToImage(mol, size=size, kekulize=True)
    bio = BytesIO()
    img.save(bio, format='PNG')
    bio.seek(0)
    return bio


def save_mol_to_sdf(mol, out_path, confId=None):
    if Chem is None:
        raise ImportError("RDKit not installed; cannot save molecule")
    writer = Chem.SDWriter(out_path)
    if confId is not None:
        writer.write(mol, confId=confId)
    else:
        writer.write(mol)
    writer.close()


def convert_to_pdbqt(input_path: str, output_path: str) -> None:
    """Convert common structure formats to PDBQT using Open Babel CLI if available.

    Raises FileNotFoundError if `obabel` not found, or subprocess.CalledProcessError on failure.
    """
    import shutil
    import subprocess
    import os

    obabel = shutil.which("obabel") or shutil.which("babel")
    if obabel:
        # Build command: obabel input -O output --partialcharge gasteiger
        cmd = [obabel, input_path, "-O", output_path, "--partialcharge", "gasteiger"]
        subprocess.run(cmd, check=True)
    else:
        # Try python Open Babel bindings (pybel) as a fallback
        try:
            # pybel exposed as `openbabel.pybel` or `pybel` in different installs
            try:
                from openbabel import pybel  # noqa: F401
                import openbabel
                pb_py = pybel
            except Exception:
                import pybel as pb_py  # type: ignore
        except Exception:
            raise FileNotFoundError(
                "Open Babel `obabel` binary not found on PATH and Open Babel Python bindings not installed; install Open Babel or use conda-forge package named openbabel."
            )

        # Infer format by extension
        fmt = os.path.splitext(input_path)[1].lstrip('.').lower()
        if fmt == 'pdbqt':
            fmt = 'pdbqt'
        elif fmt == 'pdb':
            fmt = 'pdb'
        elif fmt in ('sdf', 'mol', 'mol2'):
            # pybel recognizes 'sdf' and others
            pass
        # Read first molecule and write PDBQT
        mols = list(pb_py.readfile(fmt, input_path))
        if not mols:
            raise FileNotFoundError(f"No molecules read from {input_path} using Open Babel bindings")
        mol = mols[0]
        # Write out as pdbqt; pybel will call OBMol.Write which supports 'pdbqt' when built with OB
        mol.write('pdbqt', output_path)
    if not os.path.exists(output_path):
        raise FileNotFoundError(f"Conversion to PDBQT failed; {output_path} not created")


def pdb_to_3d_html(pdb_block: str, width: int = 400, height: int = 300) -> str:
    """Return HTML string for embedded 3D viewer using py3Dmol.
    Raises ImportError if py3Dmol not installed.
    """
    if py3Dmol is None:
        raise ImportError("py3Dmol is not available; install py3Dmol")
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_block, 'pdb')
    view.setStyle({'stick':{}})
    view.zoomTo()
    # Use internal _make_html to get full HTML
    return view._make_html()


def compute_grid_from_pdb(pdb_path: str, margin: float = 5.0):
    """Compute bounding box center and sizes from a PDB-like file.

    Returns (center_x, center_y, center_z), (size_x, size_y, size_z)
    """
    xs = []
    ys = []
    zs = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
                except Exception:
                    continue
    if not xs:
        raise ValueError('No atom coordinates found in PDB file')
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    minz, maxz = min(zs), max(zs)
    center = ((minx + maxx) / 2.0, (miny + maxy) / 2.0, (minz + maxz) / 2.0)
    size = ((maxx - minx) + margin, (maxy - miny) + margin, (maxz - minz) + margin)
    return center, size

