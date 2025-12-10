from docking.insights import compute_insights, compute_pose_interpretation, compute_contacts
import pytest


def test_compute_insights_simple():
    poses = [{'score': -10.0}, {'score': -7.0}, {'score': -5.0}]
    ins = compute_insights(poses)
    assert ins['count'] == 3
    assert ins['best'] == -10.0
    assert ins['avg'] == pytest.approx((-10.0 + -7.0 + -5.0) / 3.0)
    assert ins['range'][0] == -10.0 and ins['range'][1] == -5.0
    assert 'Strong predicted binder' in ins['interpretation']
    assert ins['promising_count'] == 2
    assert ins['strong_count'] == 1
    assert ins['promising_fraction'] == pytest.approx(2 / 3)


def test_compute_pose_interpretation_pdb_fallback():
    pdb = """
ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   UNL A   1       1.000   0.000   0.000  1.00  0.00           C
END
"""
    pose = {'pdb': pdb}
    # no receptor passed
    interp = compute_pose_interpretation(pose)
    assert 'PDB contains' in interp['text']
    assert 'atoms' in interp['details'] and interp['details']['atoms'] >= 2


def test_compute_pose_interpretation_rdkit_if_available():
    # Only run if RDKit is installed
    try:
        from rdkit import Chem
    except Exception:
        pytest.skip('RDKit not installed, skipping RDKit-based interpretation test')
    mol = Chem.MolFromSmiles('CCO')
    mol = Chem.AddHs(mol)
    pose = {'mol': mol, 'score': -8.0}
    interp = compute_pose_interpretation(pose)
    assert 'molecular_weight' in interp['details']
    assert 'ligand_efficiency' in interp['details']
    assert 'Ligand efficiency' in interp['text'] or 'LE:' in interp['text']
    # For score -8 and ethanol (3 heavy atoms), LE = -8/3 ~= -2.67 (Excellent)
    assert 'Excellent' in interp['text']


def test_compute_contacts_basic():
    receptor = """
ATOM      1  C   REC A   1       0.000   0.000   0.000  1.00  0.00           C
END
"""
    ligand = """
ATOM      1  C   LIG A   1       0.000   0.000   3.000  1.00  0.00           C
END
"""
    pose = {'pdb': ligand}
    ctr = compute_contacts(pose, receptor, cutoff=4.0)
    assert ctr['contacts'] == 1
    assert abs(ctr['closest'] - 3.0) < 1e-3
