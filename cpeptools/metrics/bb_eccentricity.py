import mdtraj as md
from rdkit import Chem
import tempfile
from ..mol_ops import get_largest_ring
from ..geometry import get_convex_hull, get_pca, get_eccentricity, fit_ellipse

def calculate_eccentricity(traj, smiles = None):
    tmp_dir = tempfile.mkdtemp()
    pdb_filename = tempfile.mktemp(suffix=".pdb", dir=tmp_dir)
    traj[0].save(pdb_filename)
    mol = Chem.MolFromPDBFile(pdb_filename, removeHs = False)

    if smiles is not None:
        ref = Chem.MolFromSmiles(smiles, sanitize = True)
        mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)

    # try: #some structures might be non-sensical and gets NaN
    indices = get_largest_ring(mol)
    for frame in traj:
        xyz = frame.xyz[0][indices, :] #only backbone indices

        xyz, variance_ratio = get_pca(xyz)
        ch_xyz = get_convex_hull(xyz)

        e_obj = fit_ellipse(xyz[:,0], xyz[:,1])
        ch_e_obj = fit_ellipse(ch_xyz[:,0], ch_xyz[:,1])
        yield get_eccentricity(e_obj), get_eccentricity(ch_e_obj), variance_ratio

    # except Exception as e: #FIXME ugly
    #     print(e)
    #     pass
