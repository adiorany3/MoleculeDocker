from docking.utils import compute_grid_from_pdb


def test_compute_grid_from_pdb():
    center, size = compute_grid_from_pdb('examples/receptor.pdb', margin=2.0)
    cx, cy, cz = center
    sx, sy, sz = size
    assert isinstance(cx, float)
    assert sx > 0
    assert sy > 0
    assert sz > 0
