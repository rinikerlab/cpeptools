from .utils import get_largest_ring_indices
from ..mol_ops import get_largest_ring
from ..geometry import get_convex_hull, get_pca, get_eccentricity, ellipse_angle_of_rotation, fit_ellipse, ellipse_center,center_2D_points, get_info_from_e_obj, place_points_on_ellipse
import numpy as np

def construct_ellipse_helper(e_obj_1, e_obj_2):
    from shapely.geometry import Polygon

    ax1,ang1,cen1 = get_info_from_e_obj(e_obj_1)
    ax2,ang2,cen2 = get_info_from_e_obj(e_obj_2)
    ellipse_1 = place_points_on_ellipse(ax1[0], ax1[1],ang1, cen1[0], cen1[1])
    ellipse_2 = place_points_on_ellipse(ax2[0], ax2[1],ang2, cen2[0], cen2[1])
    return  Polygon(np.transpose(ellipse_1)), Polygon(np.transpose(ellipse_2))

def calculate_ellipse_overlap(traj, smiles = None, metric = "Tanimoto"):

    def tanimoto(e_obj_1, e_obj_2):
        e1, e2 = construct_ellipse_helper(e_obj_1, e_obj_2)
        return e1.intersection(e2).area/e1.union(e2).area

    def norm_overlap(e_obj_1, e_obj_2):
        e1, e2 = construct_ellipse_helper(e_obj_1, e_obj_2)
        return e1.intersection(e2).area/e1.area


    metrics = {
        "Tanimoto" : tanimoto,
        "Percent" : norm_overlap, #looking at a few images, not so useful,
    }

    indices = get_largest_ring_indices(traj, smiles)

    for frame in traj:
        xyz = frame.xyz[0][indices, :] #only backbone indices

        xyz, variance_ratio = get_pca(xyz)

        e_obj = fit_ellipse(xyz[:,0], xyz[:,1])
        ch_xyz = get_convex_hull(xyz)
        ch_e_obj = fit_ellipse(ch_xyz[:,0], ch_xyz[:,1])

        yield metrics[metric](e_obj, ch_e_obj)
