import mdtraj as md
from rdkit import Chem
import tempfile
from ..mol_ops import get_largest_ring
from ..geometry import get_convex_hull, get_pca, ellipse_angle_of_rotation, fit_ellipse, ellipse_center
import numpy as np

"""
TODO:
    should merge it somehow with calculate_eccentricity
"""
def calculate_rotation(traj, smiles = None):
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

        e_obj = fit_ellipse(xyz[:,0], xyz[:,1])

        ellipse_tilt_angle =   ellipse_angle_of_rotation(e_obj) #radians, using `ellipse_angle_of_rotation2` gives large rotation angles where the actual rotation is close to none (i.e. completely almost a half turn)
        center =  ellipse_center(e_obj)
        angle_to_horizontal = np.arctan2(xyz[0,1] - center[1], xyz[0,0] - center[0]) #angle between the line of center to first point in `xyz` and the horizontal

        angle = np.absolute(angle_to_horizontal - ellipse_tilt_angle)
        # print(angle_to_horizontal - ellipse_tilt_angle)
        while angle > np.pi/2:
            angle -= np.pi/2
        yield angle

    # except Exception as e: #FIXME ugly
    #     print(e)
    #     pass
