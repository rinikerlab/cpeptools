import numpy as np
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import AllChem
from ..mol_ops import *
from ..geometry import calculate_ellipse_radii, get_points_on_ellipse
from .bmat import *
from rdkit import DistanceGeometry


def bound_matrix_from_ellipse(mol, angle, update_scheme = modify_bound_matrix, eccentricity = 0.99, bond_scale_factor = 1.0):
    """
    symmetric matrix
    """
    bmat = AllChem.GetMoleculeBoundsMatrix(mol)
    ring_indices = get_largest_ring(mol)

    bond_length_list = get_ring_bond_length_list(bmat, ring_indices, scale_factor = bond_scale_factor)
    perimeter = sum(bond_length_list)
    a, b =  calculate_ellipse_radii( guess = (perimeter/2/np.pi, perimeter/2/np.pi), eccentricity = 0.99, perimeter = perimeter) #angstroms

    coord_2d = get_points_on_ellipse(a,b, len(ring_indices), bond_length_list, startAngle = angle, verbose = False)
    bound_mat = cdist(coord_2d, coord_2d)

    bmat = update_scheme(bmat, new_matrix = bound_mat, modify = "upper", indices = ring_indices, verbose = False)
    DistanceGeometry.DoTriangleSmoothing(bmat)

    return bmat


def get_amide_pairwise_coulomb_interaction(mol, counter_charge = False, scale_factor = 0.5):
    from itertools import product
    NH = get_amine_H(mol)
    CO = get_carbonyl_O(mol)
    to_exclude = {*get_1_4_pairs(mol) , *get_1_5_pairs(mol)}

    props = AllChem.MMFFGetMoleculeProperties(mol)

    CPCIs = {}
    if counter_charge:
        CPCIs = {**CPCIs, **{(O_index, H_index) : (props.GetMMFFPartialCharge(O_index) * props.GetMMFFPartialCharge(H_index) * scale_factor if (O_index, H_index) not in to_exclude else props.GetMMFFPartialCharge(O_index) * props.GetMMFFPartialCharge(H_index) * scale_factor *-1) for O_index, H_index in product(CO, NH) }}
    else:
        CPCIs = {**CPCIs, **{(O_index, H_index) : (props.GetMMFFPartialCharge(O_index) * props.GetMMFFPartialCharge(H_index) * scale_factor) for O_index, H_index in product(CO, NH) if (O_index, H_index) not in to_exclude  }}

    # params.verbose = True
    return CPCIs
