import mdtraj as md
from scipy.spatial import Delaunay
import numpy as np
from ..geometry import *
from .utils import round_to_nearest

def in_hull(sidechain_coords, backbone_coords ):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed

    returns Indices of simplices containing each point. Points outside the triangulation get the value -1.
    """
    hull = Delaunay(backbone_coords)
    return hull.find_simplex(sidechain_coords)

def identify_interference(traj, backbone_selection = "name N or name CA or name C or name O or name H or name CN or name CB", sidechain_selection = None):
    #TODO assumes more than one frame otherwise the shape cmd would break
    #FIXME something can be neither backbone or sidechain
    idx2name = {i : val for i,val in enumerate(traj.topology.atoms)}
    coords = np.array(traj.xyz)

    if type(backbone_selection) is str:
        backbone_idx = traj.topology.select(backbone_selection)
    elif type(backbone_selection) is list and len({type(i) for i in backbone_selection}) == 1 and type(backbone_selection[0]) is int:
        backbone_idx = np.array(backbone_selection)

    if sidechain_selection is None: #all that are not backbone are considered sidechain
        if type(backbone_selection) is str:
            sidechain_idx = traj.topology.select("not ({})".format(backbone_selection))
        elif type(backbone_selection) is list and len({type(i) for i in backbone_selection}) == 1 and type(backbone_selection[0]) is int:
            sidechain_idx = np.array(list(set(range(len(idx2name))) - set(backbone_idx)))

    elif type(sidechain_selection) is str:
        sidechain_idx = traj.topology.select(sidechain_selection)
    elif type(sidechain_selection) is list and len({type(i) for i in sidechain_selection}) == 1 and type(sidechain_selection[0]) is int:
        sidechain_idx = np.array(sidechain_selection)

    backbone_coords =  coords[:,backbone_idx, :]
    sidechain_coords =  coords[:,sidechain_idx, :]
    out = {} #frame idx as key, atom idx that are in convex hull as values
    for i in range(len(traj)):
        tmp = [idx2name[j] for j in  sidechain_idx[in_hull(sidechain_coords[i], backbone_coords[i]) > 0]]
        if len(tmp) > 0:
            out[i] = tmp
    return out

#================================================

def heatmap_barchart_helper(mat, counts, lab_x, lab_y, title = ""):
    import matplotlib.pyplot as plt

    fig, ax_arr = plt.subplots(1,2,figsize=(25,15))
    ax = ax_arr[0]

    im = ax.imshow(mat,aspect="auto")

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(lab_x)))
    ax.set_yticks(np.arange(len(lab_y)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(lab_x)
    ax.set_yticklabels(lab_y)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
             rotation_mode="anchor", fontsize = 24)
    plt.setp(ax.get_yticklabels(), fontsize = 24)



    # Loop over data dimensions and create text annotations.
    for i in range(len(lab_y)):
        for j in range(len(lab_x)):
            text = ax.text(j, i, mat[i, j],
                           ha="center", va="center", color="w", size = 20)

    ax = ax_arr[1]
    ax.barh(lab_y, counts, height = 5.0, log = True, )

    plt.xticks(fontsize=24) #only changes the bar chart
    plt.gca().invert_yaxis()
    plt.grid()
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(title, fontsize=36)
    plt.show()

#FIXME
def plot_sidechain_influence_of(name, path = "/home/shuwang/Documents/Modelling/CP/Codes/"):
    with open(path + "{}.pickle".format(name), "rb") as f:
        ch = pickle.load(f)
    with open(path + "{}_eccen.pickle".format(name), "rb") as f:
        eccen = pickle.load(f)


    mat = {}
    faulty = 0
    reweighting = {e : 0 for e in range(0,105,10)}

    for i in eccen:

        for j in eccen[i]:

            eccen_value = eccen[i][j][0]
            if np.isnan(eccen_value):
                faulty+=1
                continue

            val = int(round_to_nearest(0.1)(eccen_value) * 100)
            try:
                tmp = {k.residue.name + str(k.residue.index + 1) for k in ch[i][j]}
            except KeyError:
                tmp = set(["NULL"])
            reweighting[val] += 1
            for k in tmp:
                if k in mat:
                    mat[k][val] += 1
                else:
                    mat[k] = {e : 0 for e in range(0,105,10)}
#     print(mat, faulty)
    residues = list(mat.keys())
    eccen_range = list(mat[residues[0]].keys())

    out = np.array([list(mat[i].values()) for i in mat]) #turn to matrix


    out[np.isnan(out)] = 0
    out = out.T
    out = out / (np.sum(out, 1)[:, None] + 0.0001) #turn count to conditional probability (row-wise)
    out = np.round(out , 3)
    heatmap_barchart_helper(out,  list(reweighting.values()), residues, eccen_range, name)
    del ch, eccen
