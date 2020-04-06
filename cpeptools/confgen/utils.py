import numpy as np
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import AllChem
from ..mol_ops import *
from ..geometry import calculate_ellipse_radii, get_points_on_ellipse
from .bmat import *
from rdkit import DistanceGeometry


def bound_matrix_from_ellipse(mol, angle, eccentricity = 0.99, bond_scale_factor = 1.0, update_scheme = modify_bound_matrix):
    """
    Modifies the bounds matrix entries of those atoms lying on the 
    largest ring of the input molecule. The modifications are based 
    on distances constraints derived from fitting these atoms onto 
    the circumference of an ellipse. This effectively reduces the 
    allowed conformation space spanned by the cyclic molecule as 
    defined by the bounds matrix.

    For detail see method section of https://pubs.acs.org/doi/10.1021/acs.jcim.0c00025

    Parameters
    -----------
    mol : rdkit.Chem.Mol
    angle : 
        The angle at which the first ring atom makes with the 
        horizantal. This param stipulates at which orientation the 
        ring is squashed.
    eccentricity : float
        How squashed is the ellipse that the ring atoms are fitted 
        onto. The default is 0.99. It is recommended to use a very 
        high value(eccentricity ranges from 0 to 1 for ellipse).
    bond_scale_factor : float
        The circumference of the ellipse that undergoes the fitting 
        is equal to the sum of all ring bond lengths. This param 
        controls whether the circumference is scaled. For smaller 
        ring systems, it can be good to scale up this param. The 
        default is 1.0, i.e. no scaling.
    update_scheme: function
        Rule used to modify the initial bound matrix in accordance 
        to the ellipse fit. The default function reduces both the 
        lower and upper bound whenever ellipse fit gives smaller 
        bounds.

    Return
    -------
    bmat : numpy.matrix
        The updated bound matrix, which can be subsequently fed 
        into RDKit's conformer generation workflow.
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
    """
    The heuristic used to identify amide pairwise Coulombic 
    interactions in a given molecule, i.e. those 'O and 'H' in 
    {C=O* ---- *H-N} that can in principle create hydrogen bonds, 
    which is mimicked by attractive interaction.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
    counter_charge : bool
        Whether repulsive interaction is imposed on {C=O* ---- *H-N}
        pairs that are separated by 4 or 5 atoms. Default is False, 
        which means these pairs are ignored.
    scale_factor : float
        The fraction that is multiplied to the MMFF reference value 
        for charge product. Default is 0.5, which was empirically 
        observed to work best for some cyclic peptide. It is not 
        recommended to use the full MMFF charge strength as added 
        interactions.
    
    Return
    ----------
    CPCIs : dict
        A dictionary with tuples as key and float as value, where 
        tuple enlists the atom pair to create interaction, the 
        value is the attraction strength. This dictionary can be 
        subsequently fed into RDKit's conformer generation workflow.
    """

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

    return CPCIs
