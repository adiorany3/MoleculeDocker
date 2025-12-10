from docking.utils import convert_to_pdbqt
import shutil
import pytest
import os


def test_convert_to_pdbqt_obabel_missing(tmp_path):
    # If obabel isn't installed, helper should raise FileNotFoundError
    if shutil.which('obabel') is None and shutil.which('babel') is None:
        with pytest.raises(FileNotFoundError):
            convert_to_pdbqt('examples/ligand.sdf', str(tmp_path / 'out.pdbqt'))
    else:
        # If obabel exists, we just ensure conversion runs (output file created)
        out = str(tmp_path / 'out.pdbqt')
        convert_to_pdbqt('examples/ligand.sdf', out)
        assert os.path.exists(out)
