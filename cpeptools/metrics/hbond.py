import mdtraj as md
import numpy as np
from collections import defaultdict

def calculate_hbond_occurences(traj):
    topology = traj.topology
    get_label = lambda hbond : '%s -- %s' % (topology.atom(hbond[0]), topology.atom(hbond[2]))
    # occurences = {get_label(i) : 0 for i in  md.baker_hubbard(traj, freq = 0, periodic=False)}
    occurences = defaultdict(int)

    for i in traj:
        for j in md.baker_hubbard(i, periodic = True):
            occurences[get_label(j)] += 1
        # yield occurences

    labels = [str(i) for i in topology.residues]
    return {"occurences" : dict(occurences), "labels" : labels, "num_frames" : len(traj)}


#### plotting



def plot_heatmap_helper(mat, lab_x, lab_y, title = "", save_to = None):
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(1,1,figsize=(19,15))
    # ax = ax_arr[0]

    im = ax.imshow(mat, cmap="Blues", aspect="auto")

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



    # # Loop over data dimensions and create text annotations.
    # for i in range(len(lab_y)):
    #     for j in range(len(lab_x)):
    #         text = ax.text(j, i, mat[i, j],
    #                        ha="center", va="center", color="w", size = 20)

    # fig.tight_layout()

    im.set_clim(vmin=0, vmax=1) #limiting to normalised range

    fig.suptitle(title, fontsize=36)
    plt.colorbar(im)
    plt.gcf().subplots_adjust(bottom=0.2)
    if save_to is not None:
        plt.savefig(save_to)
    else:
        plt.show()


def plot_heatmap(lookup, save_to = None, **kwargs): #calculate_hbond_occurences
    import copy
    donor_labels = copy.deepcopy(lookup["labels"])
    acceptor_labels = copy.deepcopy(lookup["labels"])
    for i in lookup["occurences"]:
        lookup["occurences"][i] /= lookup["num_frames"]  #get percentage
        donor, acceptor = i.split(" -- ")
        if donor.split("-")[-1] != "N" and donor not in donor_labels:
            donor_labels.append(donor)
        if (acceptor.split("-")[-1] != "N" and acceptor.split("-")[-1] != "O") and acceptor not in acceptor_labels:
            acceptor_labels.append(acceptor)

    mat = np.zeros((len(acceptor_labels), len(donor_labels)))

    for i in lookup["occurences"]:
        value_in_cell = lookup["occurences"][i]  #get percentage
        # print(i)
        donor, acceptor = i.split(" -- ")
        if donor.split("-")[-1] != "N" and (acceptor.split("-")[-1] != "N" and acceptor.split("-")[-1] != "O"): #both are sidechain
            mat[acceptor_labels.index(acceptor)][donor_labels.index(donor)] += value_in_cell
        elif  (acceptor.split("-")[-1] != "N" and acceptor.split("-")[-1] != "O"): #only acceptor is sidechain
            mat[acceptor_labels.index(acceptor)][donor_labels.index(donor.split("-")[0])] += value_in_cell
        elif donor.split("-")[-1] != "N": #only donor is sidechain
            mat[acceptor_labels.index(acceptor.split("-")[0])][donor_labels.index(donor)] += value_in_cell

        else: #both are NOT sidechain
            mat[acceptor_labels.index(acceptor.split("-")[0])][donor_labels.index(donor.split("-")[0])] += value_in_cell

    for idx,val in enumerate(donor_labels):
        if val.split("-")[-1] == "OG1":
            donor_labels[idx] = val.split("-")[0] + "-OH"
    for idx,val in enumerate(acceptor_labels):
        if val.split("-")[-1] == "OG1":
            acceptor_labels[idx] = val.split("-")[0] + "-OH"
    1
    plot_heatmap_helper(mat, donor_labels, acceptor_labels, kwargs["title"], save_to)
