from ..mol_ops import *
from ..geometry import get_points_on_ellipse, calculate_ellipse_radii
from .bmat import modify_bound_matrix
import copy
from scipy.spatial.distance import cdist
import numpy as np
from rdkit import Chem, DistanceGeometry
from rdkit.Chem import AllChem

def generate_conformers_with_eccentricity(mol, num_conf, angle, new_mol = True):
    print("Generating {} Conformers for {} at angle {}".format(num_conf, Chem.MolToSmiles(mol), angle))

    if new_mol:
        mol = copy.deepcopy(mol)
        mol.RemoveAllConformers()

    NH = get_amine_H(mol)
    CO = get_carbonyl_O(mol)
    to_exclude = {*get_1_4_pairs(mol) , *get_1_5_pairs(mol)}

    backbone = get_largest_ring(mol)
    perimeter = len(backbone) * 1.5 #assume bond length are ~1.5 angstrom
    a, b =  calculate_ellipse_radii( guess = (perimeter/2/np.pi, perimeter/2/np.pi), eccentricity = 0.99, perimeter = perimeter) #angstroms

    ####charging
    props = AllChem.MMFFGetMoleculeProperties(mol)

    pairCharges = {}
    # #FIXME below
    scale_factor = 0.1
    # pairCharges = {**pairCharges, **{(O_index, H_index) : (props.GetMMFFPartialCharge(O_index) * props.GetMMFFPartialCharge(H_index) * scale_factor if (O_index, H_index) not in to_exclude else props.GetMMFFPartialCharge(O_index) * props.GetMMFFPartialCharge(H_index) * scale_factor *-1) for O_index, H_index in product(CO, NH) }}
    #
    # pairCharges = {**pairCharges, **{(O_index, H_index) : (props.GetMMFFPartialCharge(O_index) * props.GetMMFFPartialCharge(H_index) * scale_factor) for O_index, H_index in product(CO, NH) if (O_index, H_index) not in to_exclude  }}
    # #FIXME above

    # params.verbose = True
    print("Pair charges : " , (pairCharges))


    amide_indices = get_atom_mapping(mol, "[O:1]=[C]@;-[NX3]-[CH3:4]") + get_atom_mapping(mol, "[O:1]=[C]@;-[NX3H1]-[H:4]") #FIXME ugly

    params = AllChem.ETKDG()
    bmat = AllChem.GetMoleculeBoundsMatrix(mol)

    coord_2d = get_points_on_ellipse(a,b, len(backbone), startAngle = angle, verbose = False)
    #FIXME below
    upper_bound_mat = cdist(coord_2d, coord_2d)
    # subset = range(0,len(backbone), 3) #TODO
    # upper_bound_mat = upper_bound_mat[[[i] for i in subset], subset]
    # backbone = [backbone[i] for i in subset]
    ##############################################
    #backbone constriants
    bmat = modify_bound_matrix(bmat,  new_matrix = upper_bound_mat, modify = "upper",  indices = backbone, verbose = False)
    #FIXME above

    DistanceGeometry.DoTriangleSmoothing(bmat)


    # trans-peptide bond constraints #FIXME
    # for x,y in amide_indices: #order matters
    #     bmat = modify_bound_matrix(bmat, value = 3.6 , modify = "lower",  indices = [x,y], check = False)
    # for x,y in amide_indices:
    #     bmat = modify_bound_matrix(bmat, value = 3.7 , modify = "upper",  indices = [x,y], check = False)


    ##FIXME below
    params.useRandomCoords = True
    # AllChem.EmbedMultipleConfs(mol, num_conf, params, pairCharges = pairCharges, boundsMatrix = bmat )
    AllChem.EmbedMultipleConfs(mol, num_conf, params, boundsMatrix = bmat)
    # # AllChem.EmbedMultipleConfs(mol, num_conf , params)
    ##FIXME above

    return mol
