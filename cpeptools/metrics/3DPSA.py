import os
import numpy as np

from typing import List, Tuple, Dict

from pymol import cmd

def do_psa3d(pdb_paths:List[str]=None, gromacs_trajs:List[Tuple[str,str]]=None, visualize:bool=True)->Dict[str, (float or List[float])]:
    """ This function allows quick 3d-PSA calculation with pymol via files.
    author: Benjamin Schroeder

    Parameters
    ----------
    pdb_paths : List[str], optional
        this argument takes a list of paths to pdb files.
    gromacs_trajs :
        this arguments takes a list of tuples containing the path to the gro file and second to the xtc file. (or pdb and dcd file, you only need to make sure topology first, traj second)
    visualize :
        visualize the results in pymol directly

    Returns
    -------
        Dict[str, (float or List[float])]
            A dictionary containing the mean, std and median + all calculated psas for all files.
    """
    # polar_surface area
    #imports
    if(visualize):
        import pymol
        pymol.finish_launching(['pymol', '-q'])   #quiet pymol '-qic'

    #loading files
    cmd.reinitialize()
    if(pdb_paths!=None):
        for pdb_path in pdb_paths:
            cmd.load(pdb_path)
    elif(gromacs_trajs!=None):
        for (gro_file, xtc_file) in gromacs_trajs:
            obj_name = os.path.splitext(os.path.basename(gro_file)[0])
            cmd.load(gro_file,object=obj_name)
            cmd.load_traj(xtc_file, object=obj_name)
    else:
        raise IOError("give me files!!!")


    ##extend Pymol with psa calc func
    cmd.extend("psa3d", _calc_psa3d)

    #calculate psa
    obj_psas = _calc_psa3d(visualization=visualize)

    return obj_psas


def _calc_psa3d(obj_list: List[str] = None,
                      default_polar_atom_elements: Tuple[str] = ("O", "P", "S", "N"), select_H_atoms_neighboring_polar_atom: bool = True,
                      upper_atom_polarity_cutoff: float = 0.5, lower_atom_polarity_cutoff: float = -0.5,
                      not_select_atom=None, additional_selection_requirement: str = "",
                      probe_sphere_radius=1.4, vdw_radii: Dict[str, float] = None,
                      state_selection: int or List[int] = -1, verbose:bool =False, visualization:bool =False) -> Dict[str, (float or List[float])]:
    """
    function to calculate the 3d polar surface area (3D-PSA) of molecules in Interface_Pymol for all the snapshots in a MD trajectory.
    author: Benjamin Schroeder

    Parameters
    ----------
    obj_list: list, optional
        list of pymol objects (Default = "cmpd1")
    default_polar_atom_elements : Tuple[float], optional
        gives the Atom elements, that are to be considered as polar atoms by default (Default = ("O", "P", "S", "N"))
    select_H_atoms_neighboring_polar_atom : bool, optional
        shall H-atoms be selected, that are neighboring a default polar atom? (Default = True)
    upper_atom_polarity_cutoff : float, optional
        if partial charges are assigned in pymol, they are compared vs this upper bound to decide if a non default polar atom should be considered polar. (Default = 0.5)
    lower_atom_polarity_cutoff : float, optional
        if partial charges are assigned in pymol, they are compared vs this lower bound to decide if a non default polar atom should be considered polar. (Default = -0.5)
    not_select_atom: str, optional
        Single atom name of the atom to remove from the selection (Default = None).
        Useful if you want to include only S or only P in the calculation of the 3D-PSA.
    additional_selection_requirement :  str, optional
        an additional requirement for the atom selection. (Default= resn LIG)
    probe_sphere_radius : float, optional
        radius of the probe sphere. (Default= 1.4)
    vdw_radii : Dict[str, float], optional
        vdw radii of the atoms in Pymol (Default=takes standard vdw radii)
    state_selection : int or List[int], optional
        select certain states/frames for 3D-PSA calculation or select all with -1 (Default= -1)

    Returns
    ----------
    psa_stats: Tuple[float, float, float]
        Values correspond to mean, standard deviation, and median of the 3D-PSA calculated over the simulation time
    """

    # IO
    if (obj_list is None):
        obj_list = cmd.get_names("objects")

    # set vdw radii of atoms
    # We do this, as the vdw radii of H is a bit old in PyMol
    if (isinstance(vdw_radii, type(None))):
        vdw_radii = vdw_std_radii = {"H": 1.09,
                                     "C": 1.7,
                                     "N": 1.55,
                                     "O": 1.52,
                                     "S": 1.8,
                                     }

    [cmd.alter("elem " + elem, "vdw=" + str(vdwr)) for elem, vdwr in vdw_radii.items()]

    if visualization:
        cmd.hide()
        cmd.show("sticks", "all")

    # probe sphere setting:
    # Allows to change the size of the probe sphere, calculating the PSA
    psa3d_set_sphere_solvent = "set solvent_radius=" + str(probe_sphere_radius) + "; set dot_density=4; set dot_solvent=1;"
    cmd.do(psa3d_set_sphere_solvent)

    # prebuild selection for relevant atoms
    ##build default polar atom selector
    polar_atoms_selector = "(elem " + "+".join(default_polar_atom_elements) + " "
    if (select_H_atoms_neighboring_polar_atom):
        polar_atoms_selector += " or (elem H and (neighbor elem " + "+".join(default_polar_atom_elements) + "))"
    polar_atoms_selector += ")"

    ##handle additional requirement
    if (isinstance(additional_selection_requirement, str) and len(additional_selection_requirement) > 0):
        additional_selection_requirement = " and " + additional_selection_requirement


    # Loop over objects
    obj_psa_dict = {}
    for obj in obj_list:
        cmd.frame(0)  # jump to first state/frame of traj

        # select states/frames, that should be used for the 3D-PSA:
        if (state_selection == -1):
            states = range(1, cmd.count_states(obj) + 1)  # get all states of the object
        elif (type(state_selection) == int):
            states = [state_selection]
        else:
            states = list(map(int, state_selection))  # list of ints is given

        ##Loop over all selected states
        psa = []
        for state in states:
            cmd.frame(state)

            ###select atoms IDs by partialCharge Filter - might be state dependent. Depends on underlying model.
            atomID_of_polarAtom = []
            cmd.iterate("(" + obj + " and (pc. < " + str(lower_atom_polarity_cutoff) + " or pc. > " + str(upper_atom_polarity_cutoff) + "))",
                        "atomID_of_polarAtom.append(ID)", space=locals())
            partialCSelection = " or ID " + "+".join(map(str, atomID_of_polarAtom)) if (len(atomID_of_polarAtom) > 0) else ""

            ###select all needed atoms by partialCSelection, default polar, etc ...
            if not_select_atom != None and isinstance(not_select_atom, str):
                select_string = obj + " and ((" + polar_atoms_selector + " " + additional_selection_requirement + " and not name {}".format(
                    not_select_atom) + " ) " + partialCSelection + ")"
            else:
                select_string = obj + " and " + "((" + polar_atoms_selector + " " + additional_selection_requirement + " ) " + partialCSelection + ")"

            if(verbose): print("Selection: ", select_string)
            cmd.select("polarAtoms", select_string)
            ###calc surface area
            psa.append(float(cmd.get_area("polarAtoms", state=state)))

            ###visualise Selection
            if visualization:  # visualises what the script does in RestraintMaker_Pymol
                vis = "color grey, " + obj + "\n show surface, " + obj + " and state " + str(
                    state) + "\n set transparency, 0.4, " + obj + "\n show spheres, polarAtoms\n color red, polarAtoms\n"
                cmd.do(vis)
        if(verbose): print("All PSAs : ", psa)

        ###gather data
        obj_psa_dict.update({obj: {"avg_psa": np.mean(psa), "std_psa": np.std(psa), "median_psa": np.median(psa),"all_psas":psa}})

    return obj_psa_dict



if __name__ == "__main__":
    import pymol
    pymol.finish_launching(['pymol', '-q'])  # quiet pymol '-qic'

    tmp_pdb = "./1VPR.pdb"  #not the best molecule, but well
    cmd.fetch("1VPR")
    cmd.save(tmp_pdb)

    psas = do_psa3d(pdb_paths=[tmp_pdb], visualize=True)
    print("Result: ", psas)

    os.remove(tmp_pdb)
