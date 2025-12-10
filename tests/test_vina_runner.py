from docking.vina_runner import parse_vina_output
import tempfile

VINA_EXAMPLE = """
REMARK VINA RESULT:     -6.8      0.000      0.000
MODEL        1
ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00  0.00      A
ENDMDL
REMARK VINA RESULT:     -5.2      0.000      0.000
MODEL        2
ATOM      1  C   UNL A   1       1.000   1.000   1.000  1.00  0.00      A
ENDMDL
"""


def test_parse_vina_output_tempfile():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False) as tmp:
        tmp.write(VINA_EXAMPLE)
        tmp.flush()
        path = tmp.name
    poses = parse_vina_output(path)
    assert len(poses) == 2
    assert poses[0]['score'] == -6.8
    assert 'pdbqt' in poses[0]
    assert poses[1]['score'] == -5.2

