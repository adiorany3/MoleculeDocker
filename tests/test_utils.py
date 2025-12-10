from docking.utils import allowed_file_ext, pdb_to_3d_html
import pytest


def test_allowed_file_ext():
    assert allowed_file_ext('model.pdb')
    assert allowed_file_ext('ligand.sdf')
    assert not allowed_file_ext('image.png')


def test_pdb_to_3d_html_skip_if_missing():
    # Import or skip if py3Dmol not installed
    try:
        import py3Dmol
    except Exception:
        pytest.skip('py3Dmol not installed, skipping 3D viewer test')

    pdb_block = """
ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00  0.00           C
END
"""
    html = pdb_to_3d_html(pdb_block)
    assert '<html' in html.lower()
