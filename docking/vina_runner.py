import shutil
import subprocess
import re
import tempfile
import os

VINA_EXEC_NAME = "vina"


def is_vina_available():
    """Return True if `vina` is on PATH."""
    return shutil.which(VINA_EXEC_NAME) is not None


def run_vina(receptor_pdbqt, ligand_pdbqt, center, size, out="out.pdbqt", log="vina.log", exhaustiveness=8, num_modes=9):
    """Run AutoDock Vina. Inputs must be PDBQT.

    Parameters
    - receptor_pdbqt, ligand_pdbqt: paths to PDBQT files
    - center: (x,y,z)
    - size: (x,y,z)
    - out: output path for poses
    - log: log file

    Raises subprocess.CalledProcessError on failure.
    """
    if not is_vina_available():
        raise FileNotFoundError("Vina executable not found on PATH")

    cmd = [
        VINA_EXEC_NAME,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--out", out,
        "--log", log,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
    ]

    subprocess.run(cmd, check=True)


def parse_vina_output(pdbqt_path):
    """Parse PDBQT output file into a list of poses (score & pdbqt string).

    Returns list of dicts: {'score': float, 'pdbqt': str}
    """
    poses = []
    if not os.path.exists(pdbqt_path):
        return poses

    with open(pdbqt_path, "r") as f:
        lines = f.readlines()

    current_pose = []
    score = None
    for line in lines:
        if line.startswith("REMARK VINA RESULT:"):
            # Example: REMARK VINA RESULT:     -7.2      0.000      0.000
            m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", line)
            if m:
                score = float(m.group(1))
        if line.startswith("MODEL"):
            current_pose = [line]
        elif line.startswith("ENDMDL"):
            current_pose.append(line)
            poses.append({"score": score, "pdbqt": "".join(current_pose)})
            current_pose = []
            score = None
        else:
            if current_pose is not None:
                current_pose.append(line)

    return poses
