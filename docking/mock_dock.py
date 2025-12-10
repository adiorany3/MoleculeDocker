try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception as e:
    Chem = None
    AllChem = None

import os
import random


def read_molecule(path):
    """Try common file formats and return an RDKit Mol if possible, else None."""
    if Chem is None:
        raise ImportError("RDKit not installed; please install rdkit-pypi or use conda to install rdkit")

    ext = os.path.splitext(path)[1].lower()
    if ext in [".sdf", ".mol"]:
        mol = Chem.MolFromMolFile(path, removeHs=False)
    elif ext in [".pdb", ".pdbqt"]:
        mol = Chem.MolFromPDBFile(path, removeHs=False)
    else:
        mol = None
    return mol


def mock_dock(ligand_path, num_poses=5):
    """Perform a mock docking by generating conformers and scoring them with UFF energies.

    Returns list of dicts: {'score': float, 'mol': RDKit Mol}
    """
    # If RDKit is available, generate conformers and return RDKit Mol-based poses
    if Chem is not None and AllChem is not None:
        mol = read_molecule(ligand_path)
        if mol is None:
            raise ValueError("Could not read ligand file. Use SDF, MOL, or PDB formats.")

        mol = Chem.AddHs(mol)
        # Generate conformers
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_poses, randomSeed=42)
        # Optimize conformers
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=200)
        scored = []
        for confId, res in zip(cids, results):
            # res is a tuple (status, energy)
            if isinstance(res, tuple):
                energy = res[1]
            else:
                energy = float(res)
            scored.append({"score": float(energy), "mol": mol, "confId": int(confId)})

        # sort by lowest energy
        scored.sort(key=lambda x: x["score"]) 
        # keep top num_poses
        return scored[:num_poses]

    # If RDKit is not available, return simple text-based PDB mock poses (no images)
    import random
    poses = []
    for i in range(num_poses):
        score = float(random.uniform(-10.0, 5.0))
        # Build a minimal PDB block with 6 carbon atoms placed randomly
        atoms = []
        for a in range(1, 7):
            x = random.uniform(-2.0, 2.0)
            y = random.uniform(-2.0, 2.0)
            z = random.uniform(-2.0, 2.0)
            atoms.append(f"HETATM{a:5d}  C   UNL A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
        pdb_block = "\n".join(atoms) + "\nEND\n"
        poses.append({"score": score, "pdb": pdb_block})
    return poses
