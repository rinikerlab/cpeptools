import mdtraj as md
from scipy.spatial import Delaunay
import numpy as np
from ..geometry import *
#================================================
def in_hull(sidechain_coords, backbone_coords ):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    hull = Delaunay(backbone_coords)
    return hull.find_simplex(sidechain_coords)

def identify_interference(traj):
    #TODO assumes more than one frame otherwise the shape cmd would break
    selection = "name N or name CA or name C or name O or name H or name CN or name CB"
    idx2name = {i : val for i,val in enumerate(traj.topology.atoms)}
    coords = np.array(traj.xyz)
    backbone_idx = traj.topology.select(selection)
    sidechain_idx = traj.topology.select("not ({})".format(selection))


    backbone_coords =  np.squeeze(coords[:,[backbone_idx], :], 1)
    sidechain_coords =  np.squeeze(coords[:,[sidechain_idx], :], 1)
    out = {}
    print(backbone_coords.shape)
    for i in range(len(traj)):
        tmp = [idx2name[j] for j in  sidechain_idx[in_hull(sidechain_coords[i], backbone_coords[i]) > 0]]
        if len(tmp) > 0:
            out[i] = tmp
    return out
