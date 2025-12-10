from statistics import mean, stdev

try:
    from rdkit.Chem import Descriptors, Lipinski
    from rdkit.Chem import rdMolDescriptors
    RDKit_available = True
except Exception:
    RDKit_available = False


def compute_insights(poses: list) -> dict:
    scores = [p.get('score') for p in poses if p.get('score') is not None]
    insights = {"count": len(poses), "best": None, "avg": None, "stdev": None, "range": None, "interpretation": '', "promising_count": 0, "strong_count": 0, "promising_fraction": 0.0}
    if not scores:
        return insights
    best = min(scores)
    avg = mean(scores)
    sdev = stdev(scores) if len(scores) > 1 else 0.0
    rmin, rmax = min(scores), max(scores)
    insights.update({"best": best, "avg": avg, "stdev": sdev, "range": (rmin, rmax)})
    insights['promising_count'] = sum(1 for s in scores if s <= -7)
    insights['strong_count'] = sum(1 for s in scores if s <= -9)
    insights['promising_fraction'] = (insights['promising_count'] / len(scores)) if scores else 0.0
    # Simple heuristic
    if best <= -9:
        insights['interpretation'] = 'Strong predicted binder — consider experimental validation (🔬)'
    elif best <= -7:
        insights['interpretation'] = 'Moderate predicted binder — interesting candidate (⭐)'
    elif best <= -5:
        insights['interpretation'] = 'Weak predicted binder — may need optimization (⚠️)'
    else:
        insights['interpretation'] = 'Poor predicted binding — likely not promising (❌)'
    return insights


def compute_pose_interpretation(pose: dict, receptor_pdb_text: str = None) -> dict:
    interp = {"text": "No interpretation available", "details": {}}
    # RDKit mol descriptors
    try:
        if RDKit_available and 'mol' in pose and pose['mol'] is not None:
            mol = pose['mol']
            score = pose.get('score')
            mw = float(Descriptors.MolWt(mol))
            heavy = int(rdMolDescriptors.CalcNumHeavyAtoms(mol))
            hbd = int(Lipinski.NumHDonors(mol))
            hba = int(Lipinski.NumHAcceptors(mol))
            tpsa = float(rdMolDescriptors.CalcTPSA(mol))
            rot = int(Lipinski.NumRotatableBonds(mol))
            le = None
            if score is not None and heavy > 0:
                le = float(score) / float(heavy)

            details = {
                'molecular_weight': mw,
                'heavy_atoms': heavy,
                'h_bond_donors': hbd,
                'h_bond_acceptors': hba,
                'tpsa': tpsa,
                'rotatable_bonds': rot,
                'ligand_efficiency': le,
            }
            lines = []
            lines.append(f"MW: {mw:.1f} g/mol; Heavy atoms: {heavy}")
            lines.append(f"HBD: {hbd}, HBA: {hba}, TPSA: {tpsa:.1f}")
            lines.append(f"Rotatable bonds: {rot}")
            if le is not None:
                lines.append(f"Ligand efficiency (kcal/mol per heavy atom): {le:.3f}")
                if le <= -0.4:
                    lines.append("LE: Excellent (<= -0.4) — efficient binder per heavy atom 🔥")
                elif le <= -0.3:
                    lines.append("LE: Good (<= -0.3) — promising lead candidate ⭐")
                else:
                    lines.append("LE: Low — may need optimization ⚠️")
            if tpsa > 140:
                lines.append("TPSA high (>140): may indicate low permeability 🚫")
            interp['text'] = ' | '.join(lines)
            interp['details'] = details
            if receptor_pdb_text:
                try:
                    ctr = compute_contacts(pose, receptor_pdb_text)
                    interp['details'].update({'contacts_to_receptor': ctr})
                    if ctr.get('contacts') is not None and ctr.get('closest') is not None:
                        interp['text'] += f" | Contacts: {ctr.get('contacts')} (closest {ctr.get('closest'):.2f} Å)"
                except Exception:
                    pass
            return interp
    except Exception:
        pass

    # Fallback: PDB/PDBQT text heuristics
    try:
        pdb_text = pose.get('pdb') or pose.get('pdbqt')
        if pdb_text:
            xs = []
            ys = []
            zs = []
            atoms = 0
            for line in pdb_text.splitlines():
                if line.startswith(('ATOM', 'HETATM')):
                    atoms += 1
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        xs.append(x); ys.append(y); zs.append(z)
                    except Exception:
                        continue
            if atoms:
                span = max(xs) - min(xs) if xs else 0.0
                interp['text'] = f"PDB contains {atoms} atoms; size ≈ {span:.2f} Å"
                interp['details'] = {'atoms': atoms, 'span': span}
                if receptor_pdb_text:
                    try:
                        ctr = compute_contacts(pose, receptor_pdb_text)
                        interp['details'].update({'contacts_to_receptor': ctr})
                        if ctr.get('contacts') is not None and ctr.get('closest') is not None:
                            interp['text'] += f" | Contacts: {ctr.get('contacts')} (closest {ctr.get('closest'):.2f} Å)"
                    except Exception:
                        pass
                return interp
    except Exception:
        pass

    return interp


def compute_contacts(pose: dict, receptor_pdb_text: str, cutoff: float = 4.0) -> dict:
    """Compute simple contact metrics between receptor and ligand pose.

    Returns dict with {'contacts': int, 'closest': float, 'mean_min_distance': float}.
    """
    import numpy as _np

    # Parse receptor coordinates
    rx = []
    ry = []
    rz = []
    for line in receptor_pdb_text.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                rx.append(x); ry.append(y); rz.append(z)
            except Exception:
                continue
    if not rx:
        return {'contacts': 0, 'closest': None, 'mean_min_distance': None}
    R = _np.array([rx, ry, rz]).T

    # Parse ligand coordinates
    lx = []
    ly = []
    lz = []
    if pose.get('pdb') or pose.get('pdbqt'):
        pdb_text = pose.get('pdb') or pose.get('pdbqt')
        for line in pdb_text.splitlines():
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    lx.append(x); ly.append(y); lz.append(z)
                except Exception:
                    continue
    elif pose.get('mol') is not None:
        try:
            # RDKit mol with 3D conformation — use first conformer or confId
            from rdkit.Chem import rdMolTransforms
            from rdkit import Chem
            mol = pose['mol']
            confId = pose.get('confId')
            conf = mol.GetConformer(confId) if confId is not None else mol.GetConformer(0)
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                lx.append(pos.x); ly.append(pos.y); lz.append(pos.z)
        except Exception:
            pass
    if not lx:
        return {'contacts': 0, 'closest': None, 'mean_min_distance': None}
    L = _np.array([lx, ly, lz]).T

    # compute distance matrix
    dmat = _np.sqrt(((L[:, None, :] - R[None, :, :]) ** 2).sum(axis=2))
    # minimal distances from each ligand atom to receptor
    min_dist = dmat.min(axis=1)
    contacts = int((min_dist <= cutoff).sum())
    closest = float(min_dist.min())
    mean_min = float(min_dist.mean())
    return {'contacts': contacts, 'closest': closest, 'mean_min_distance': mean_min}
