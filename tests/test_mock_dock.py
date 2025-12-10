from docking.mock_dock import mock_dock


def test_mock_dock_sdf():
    poses = mock_dock('examples/ligand.sdf', num_poses=2)
    assert len(poses) == 2
    for p in poses:
        assert 'score' in p
        # Either RDKit mol is present or pdb text is present
        assert ('mol' in p and p['mol'] is not None) or ('pdb' in p and p['pdb'] is not None)
