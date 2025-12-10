import streamlit as st
import tempfile
import os
import streamlit.components.v1 as components
from docking.vina_runner import is_vina_available, run_vina, parse_vina_output
from docking.mock_dock import mock_dock
from docking.utils import load_molecule_image, save_mol_to_sdf, allowed_file_ext, convert_to_pdbqt, pdb_to_3d_html, compute_grid_from_pdb
try:
    import rdkit
    RDKit_available = True
except Exception:
    RDKit_available = False

from docking.insights import compute_insights, compute_pose_interpretation

st.set_page_config(page_title="Molecule Docking Online", layout="wide")

st.title("Molecule Docking Online")

st.sidebar.header("Inputs")
receptor_file = st.sidebar.file_uploader("Upload receptor (PDB or PDBQT)", type=["pdb", "pdbqt"])
ligand_file = st.sidebar.file_uploader("Upload ligand (SDF, MOL, PDB, PDBQT)", type=["sdf", "mol", "pdb", "pdbqt"]) 

# Cache uploaded receptor bytes in session_state for callback use
if receptor_file is not None:
    rb = receptor_file.getvalue()
    if st.session_state.get('receptor_bytes') != rb:
        st.session_state['receptor_bytes'] = rb
        st.session_state['receptor_filename'] = receptor_file.name

use_vina = is_vina_available()

st.sidebar.markdown(f"**Vina available:** {'✅ Yes — AutoDock Vina backend' if use_vina else '⚠️ No — mock docking fallback'}")
st.sidebar.markdown(f"**RDKit available:** {'✅ Yes — molecule rendering & SDF' if RDKit_available else '❌ No — images may not render'}")


## insights functions previously in app.py have been moved to docking/insights.py

# compute_pose_interpretation now lives in docking/insights.py and is imported at top of this file

# Compute grid button is intentionally placed before the number inputs so the handler sets session_state
# before the number_input widgets are instantiated (avoids session_state modification error).
def compute_grid_callback():
    # read bytes from session state
    if 'receptor_bytes' not in st.session_state:
        st.session_state['compute_grid_error'] = 'Please upload a receptor file first.'
        return
    import tempfile, os
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = st.session_state.get('receptor_filename', 'receptor.pdb')
            rec_path = os.path.join(tmpdir, fname)
            with open(rec_path, 'wb') as f:
                f.write(st.session_state['receptor_bytes'])
            center, size = compute_grid_from_pdb(rec_path)
            cx, cy, cz = center
            sx, sy, sz = size
            st.session_state['center_x'] = float(cx)
            st.session_state['center_y'] = float(cy)
            st.session_state['center_z'] = float(cz)
            st.session_state['size_x'] = float(sx)
            st.session_state['size_y'] = float(sy)
            st.session_state['size_z'] = float(sz)
            st.session_state['compute_grid_success'] = True
    except Exception as exc:
        st.session_state['compute_grid_error'] = str(exc)
        return

compute_grid_button = st.sidebar.button("Compute grid from receptor", on_click=compute_grid_callback)

center_x = st.sidebar.number_input("Center X", value=st.session_state.get("center_x", 0.0), step=0.1, key="center_x")
center_y = st.sidebar.number_input("Center Y", value=st.session_state.get("center_y", 0.0), step=0.1, key="center_y")
center_z = st.sidebar.number_input("Center Z", value=st.session_state.get("center_z", 0.0), step=0.1, key="center_z")
size_x = st.sidebar.number_input("Size X (Å)", value=st.session_state.get("size_x", 20.0), step=1.0, key="size_x")
size_y = st.sidebar.number_input("Size Y (Å)", value=st.session_state.get("size_y", 20.0), step=1.0, key="size_y")
size_z = st.sidebar.number_input("Size Z (Å)", value=st.session_state.get("size_z", 20.0), step=1.0, key="size_z")

num_poses = st.sidebar.slider("Number of poses", min_value=1, max_value=20, value=5)
exhaustiveness = st.sidebar.slider("Exhaustiveness (Vina)", min_value=1, max_value=32, value=8)
show_3d = st.sidebar.checkbox("Show 3D viewer (py3Dmol)", value=False)
image_width = st.sidebar.number_input("Image width (px)", min_value=100, max_value=1200, value=300, step=50, help="Width in pixels used for pose images; defaults to 300")

run_button = st.sidebar.button("Run Docking")

convert_button = st.sidebar.button("Convert to PDBQT (Open Babel)")

if run_button:
    if not receptor_file or not ligand_file:
        st.error("Please upload both receptor and ligand files.")
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save uploaded files to temp dir
            rec_path = os.path.join(tmpdir, receptor_file.name)
            lig_path = os.path.join(tmpdir, ligand_file.name)
            with open(rec_path, "wb") as f:
                f.write(receptor_file.getbuffer())
            with open(lig_path, "wb") as f:
                f.write(ligand_file.getbuffer())

            # If Vina is available, ensure receptor/ligand are PDBQT. Convert automatically if possible.
            if use_vina:
                rec_to_use = rec_path
                lig_to_use = lig_path
                # convert if needed
                if not receptor_file.name.lower().endswith('.pdbqt'):
                    try:
                        rec_to_use = os.path.join(tmpdir, os.path.splitext(receptor_file.name)[0] + '.pdbqt')
                        convert_to_pdbqt(rec_path, rec_to_use)
                        st.info('Receptor converted to PDBQT ✅')
                    except Exception as e:
                        st.warning(f'Receptor conversion failed: {e}; running mock docking')
                        use_vina = False
                if not ligand_file.name.lower().endswith('.pdbqt'):
                    try:
                        lig_to_use = os.path.join(tmpdir, os.path.splitext(ligand_file.name)[0] + '.pdbqt')
                        convert_to_pdbqt(lig_path, lig_to_use)
                        st.info('Ligand converted to PDBQT ✅')
                    except Exception as e:
                        st.warning(f'Ligand conversion failed: {e}; running mock docking')
                        use_vina = False
            if use_vina and (receptor_file.name.lower().endswith(".pdbqt") and ligand_file.name.lower().endswith(".pdbqt")):
                st.info("Using AutoDock Vina backend")
                out_pdbqt = os.path.join(tmpdir, "out.pdbqt")
                log = os.path.join(tmpdir, "vina.log")
                try:
                    run_vina(rec_to_use, lig_to_use, (center_x, center_y, center_z), (size_x, size_y, size_z), out=out_pdbqt, log=log, exhaustiveness=exhaustiveness, num_modes=num_poses)
                    poses = parse_vina_output(out_pdbqt)
                except Exception as e:
                    st.error(f"Vina failed: {e}")
                    poses = []
            else:
                st.info("Using mock docking backend (RDKit-based) — this is a demo fallback")
                # convert to SDF or use RDKit directly
                try:
                    poses = mock_dock(lig_path, num_poses)
                except Exception as e:
                    st.error(f"Mock docking failed: {e}")
                    poses = []

            if not poses:
                st.warning("No poses found")
            else:
                st.success(f"Found {len(poses)} poses")
                # Compute and display insights summary
                insights = compute_insights(poses)
                if insights.get('best') is not None:
                    best = insights['best']
                    avg = insights['avg']
                    sdev = insights['stdev']
                    rmin, rmax = insights['range']
                    interpretation = insights['interpretation']
                    st.markdown("---")
                    st.markdown("**Summary** 📊")
                    st.write(f"- Count: **{insights['count']}**")
                    # show metrics side-by-side
                    m1, m2, m3 = st.columns(3)
                    m1.metric("Best (kcal/mol)", f"{best:.3f}")
                    m2.metric("Average (kcal/mol)", f"{avg:.3f}")
                    m3.metric("Std dev", f"{sdev:.3f}" if sdev is not None else "0.000")
                    st.write(f"- Score range: **{rmin:.3f} — {rmax:.3f} kcal/mol**")
                    if sdev is not None:
                        pass
                    st.info(f"Interpretation: {interpretation}")
                    st.caption("Note: Vina-style docking scores are in kcal/mol; lower (more negative) indicates stronger predicted binding 🧭")
                    # Suggestions to improve results
                    suggestions = []
                    if insights.get('promising_fraction', 0.0) < 0.25:
                        suggestions.append('Increase `exhaustiveness` in the sidebar to allow more thorough sampling (e.g. 16-32).')
                        suggestions.append('Increase `Number of poses` to sample more conformers (e.g. 10-50).')
                        suggestions.append('Ensure receptor is prepared (PDBQT) and grid box centers on active site: use `Compute grid from receptor`.')
                    if insights.get('stdev', 0.0) > 1.5:
                        suggestions.append('High score variability: consider clustering poses or using consensus scoring for robust selection.')
                    # ligand-based suggestions
                    best_le = None
                    if 'best' in insights and insights['best'] is not None:
                        # If we have best pose and RDKit mol, compute LE in details
                        if 'details' in locals():
                            pass
                    if suggestions:
                        st.markdown("**Suggestions to improve results**")
                        for s in suggestions:
                            st.write(f"- {s}")
                cols = st.columns(3)
                # mark index of top pose (top = lowest score)
                best_score = None
                for p in poses:
                    if p.get('score') is not None:
                        if best_score is None or p['score'] < best_score:
                            best_score = p['score']
                # Find the best pose to interpret
                best_pose = None
                if best_score is not None:
                    for p in poses:
                        if p.get('score') is not None and abs(p['score'] - best_score) < 1e-6:
                            best_pose = p
                            break
                # If best pose found, compute interpretation and display
                if best_pose is not None:
                    receptor_pdb_text = None
                    if st.session_state.get('receptor_bytes'):
                        try:
                            receptor_pdb_text = st.session_state.get('receptor_bytes').decode('utf-8')
                        except Exception:
                            try:
                                receptor_pdb_text = st.session_state.get('receptor_bytes').decode('latin-1')
                            except Exception:
                                receptor_pdb_text = None
                    interp = compute_pose_interpretation(best_pose, receptor_pdb_text)
                    st.markdown("---")
                    st.markdown("**Best pose interpretation** 🔍")
                    with st.expander("Show interpretation details", expanded=False):
                        st.write(interp.get('text'))
                        details = interp.get('details') or {}
                        if details:
                            for k, v in details.items():
                                st.write(f"- **{k}**: {v}")
                    # Recommendations based on top pose details
                    more_suggestions = []
                    det = interp.get('details') or {}
                    le = None
                    if det.get('ligand_efficiency') is not None:
                        try:
                            le = float(det.get('ligand_efficiency'))
                        except Exception:
                            le = None
                    if le is not None:
                        if le <= -0.4:
                            more_suggestions.append('Excellent ligand efficiency — consider experimental validation and prioritize this compound. ✔️')
                        elif le <= -0.3:
                            more_suggestions.append('Good ligand efficiency — select for further optimization and SAR exploration. ⭐')
                        else:
                            more_suggestions.append('Low ligand efficiency — consider small-group optimization, reduce polar surface area, or improve interactions. ⚠️')
                    # contacts
                    if det.get('contacts_to_receptor'):
                        try:
                            contacts = det['contacts_to_receptor']['contacts']
                            if contacts < 3:
                                more_suggestions.append('Few close contacts with receptor detected — check grid box and ensure ligand is placed in pocket. Consider increasing `exhaustiveness` or `num_poses`.')
                        except Exception:
                            pass
                    if more_suggestions:
                        st.markdown('**Top-pose suggestions**')
                        for s in more_suggestions:
                            st.write(f"- {s}")
                for i, pose in enumerate(poses, 1):
                    col = cols[(i - 1) % 3]
                    # pose should be a dict: {'score': float, 'mol': rdkit Mol or pdbqt path}
                    score = pose.get("score")
                    mol = pose.get("mol")
                    confId = pose.get("confId")
                    pdbqt_block = pose.get("pdbqt")
                    pdb_block = pose.get("pdb")
                    img = None
                    if mol is not None:
                        try:
                            img = load_molecule_image(mol, confId=confId)
                        except Exception:
                            img = None
                    with col:
                        title = f"Pose {i}"
                        badge = ''
                        if pose.get('score') is not None:
                            s = pose.get('score')
                            if s <= -9:
                                badge = ' 🔥'
                            elif s <= -7:
                                badge = ' ⭐'
                            elif s <= -5:
                                badge = ' ⚠️'
                            else:
                                badge = ' ❌'
                        if pose.get('score') is not None and best_score is not None and abs(pose.get('score') - best_score) < 1e-6:
                            title += " — 🏆 Best"
                        title += badge
                        st.subheader(title)
                        st.write(f"Score: {score:.3f} kcal/mol" if score is not None else "Score: N/A")
                        if img is not None:
                            st.image(img, width=image_width)
                        elif pdbqt_block is not None:
                            st.code('\n'.join(pdbqt_block.splitlines()[:12]) + '\n...')
                            if show_3d:
                                try:
                                    html = pdb_to_3d_html(pdbqt_block)
                                    components.html(html, height=320)
                                except Exception:
                                    pass
                        elif pdb_block is not None:
                            st.code('\n'.join(pdb_block.splitlines()[:12]) + '\n...')
                            if show_3d:
                                try:
                                    html = pdb_to_3d_html(pdb_block)
                                    components.html(html, height=320)
                                except Exception:
                                    pass
                        # write downloadable file
                        # Save pose: prefer SDF if RDKit mol present, otherwise save PDBQT
                        if mol is not None:
                            out_path = os.path.join(tmpdir, f"pose_{i}.sdf")
                            try:
                                save_mol_to_sdf(mol, out_path, confId=confId)
                                with open(out_path, "rb") as f:
                                    st.download_button(f"⬇️ Download pose {i} (SDF)", data=f, file_name=f"pose_{i}.sdf")
                            except Exception:
                                pass
                        elif pdbqt_block is not None:
                            out_path = os.path.join(tmpdir, f"pose_{i}.pdbqt")
                            try:
                                with open(out_path, "w") as f:
                                    f.write(pdbqt_block)
                                with open(out_path, "rb") as f:
                                    st.download_button(f"⬇️ Download pose {i} (PDBQT)", data=f, file_name=f"pose_{i}.pdbqt")
                            except Exception:
                                pass
                        elif pdb_block is not None:
                            out_path = os.path.join(tmpdir, f"pose_{i}.pdb")
                            try:
                                with open(out_path, "w") as f:
                                    f.write(pdb_block)
                                with open(out_path, "rb") as f:
                                    st.download_button(f"⬇️ Download pose {i} (PDB)", data=f, file_name=f"pose_{i}.pdb")
                            except Exception:
                                pass
                        


if convert_button:
    if not receptor_file or not ligand_file:
        st.error("Please upload both receptor and ligand files to convert.")
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            rec_path = os.path.join(tmpdir, receptor_file.name)
            lig_path = os.path.join(tmpdir, ligand_file.name)
            with open(rec_path, "wb") as f:
                f.write(receptor_file.getbuffer())
            with open(lig_path, "wb") as f:
                f.write(ligand_file.getbuffer())
            try:
                rec_out = os.path.join(tmpdir, os.path.splitext(receptor_file.name)[0] + ".pdbqt")
                lig_out = os.path.join(tmpdir, os.path.splitext(ligand_file.name)[0] + ".pdbqt")
                convert_to_pdbqt(rec_path, rec_out)
                convert_to_pdbqt(lig_path, lig_out)
                st.success("Conversion completed ✅ — files are available for download")
                with open(rec_out, "rb") as f:
                    st.download_button("⬇️ Download receptor (PDBQT)", data=f, file_name=os.path.basename(rec_out))
                with open(lig_out, "rb") as f:
                    st.download_button("⬇️ Download ligand (PDBQT)", data=f, file_name=os.path.basename(lig_out))
            except FileNotFoundError as e:
                st.error(str(e))
                st.info("Install Open Babel (obabel) with:\n- conda: `conda install -c conda-forge openbabel`\n- macOS (Homebrew): `brew install open-babel`\n- Debian/Ubuntu: `sudo apt-get install openbabel`")
            except Exception as e:
                st.error(f"Conversion failed: {e}")

if st.session_state.get('compute_grid_error'):
    st.sidebar.error(st.session_state.pop('compute_grid_error'))
if st.session_state.get('compute_grid_success'):
    st.sidebar.success('Grid computed from receptor and applied — refresh done')

st.sidebar.markdown("---")
st.sidebar.markdown("Made with ❤️ — By Galuh Adi Insani")
