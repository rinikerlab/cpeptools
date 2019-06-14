#FIXME
import mdtraj as md
import numpy as np
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
from ..mol_ops import get_largest_ring
# def get_frames( which_clusters, rmsd_threshold = 0.1):
#     output = []
#     num_traj = len(traj_list)
#     for i in zip(traj_list, range(num_traj)):
#         traj = i[0]
#         indices = [int(j) -1  for j in np.load("./Trace_{}.npy".format(i[1]))]
#
#         output.append(traj[[j in which_clusters for j in indices]])
#         if i[1] == num_traj - 1:
#             # print("here", [len(k) for k in output])
#             traj = reduce(lambda a,b : a+b, output)
#             return traj.superpose(traj, 0 , atom_indices = traj.topology.select("name CA or name C or name N or name O"))


def _get_largest_ring_indices(traj, smiles = None):
    tmp_dir = tempfile.mkdtemp()
    pdb_filename = tempfile.mktemp(suffix=".pdb", dir=tmp_dir)
    traj[0].save(pdb_filename)

    if smiles is not None:
        mol = Chem.MolFromPDBFile(pdb_filename, removeHs = True)
        ref = Chem.MolFromSmiles(smiles, sanitize = True)
        mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)
    else:
        mol = Chem.MolFromPDBFile(pdb_filename, removeHs = False)

    # try: #some structures might be non-sensical and gets NaN
    indices = get_largest_ring(mol)
    
    return indices

def remove_similar_frames(traj, rmsd_threshold = 0.1):
    counter = 0
    while True:
        if len(traj) <= counter:
            break
        rmsd = md.rmsd(traj, traj, counter, atom_indices = traj.topology.select("name CA or name C or name N or name O") )
        traj = traj[(rmsd > rmsd_threshold) | (rmsd == 0)]

        counter += 1
    counter = -1
    while True:
        if len(traj) <= abs(counter):
            break
        rmsd = md.rmsd(traj, traj, len(traj) + counter )
        traj = traj[(rmsd > rmsd_threshold) | (rmsd == 0)]
        counter -= 1

    # reorder the traj based on sum of all pair rmsd
    rank = [-sum(md.rmsd(traj, traj, i)) for i in range(len(traj))]
    traj = traj[np.argsort(rank)]

    print(len(traj), " conformers to write out")
    return traj

def random_select_frames(traj, num_frames = 5):
    return traj[np.random.randint(len(traj), size = num_frames)]
