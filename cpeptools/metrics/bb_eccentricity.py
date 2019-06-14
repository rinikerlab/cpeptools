import numpy as np
from .utils import _get_largest_ring_indices
from ..geometry import get_convex_hull, get_pca, get_eccentricity, ellipse_angle_of_rotation, fit_ellipse, ellipse_center,center_2D_points, get_info_from_e_obj


def get_angle_helper(e_obj, xyz):
    """
        Angles are in range of (0,180) 0-pi

        Using convention of the projected ellipse with points ordered in **ANTICLOCKWISE** fashion

        to flip ellipse from clockwise to anticlockwise: minus 180 degrees
    """
    ellipse_tilt_angle =   ellipse_angle_of_rotation(e_obj) #radians, using `ellipse_angle_of_rotation2` gives large rotation angles where the actual rotation is close to none (i.e. completed almost a half turn)
    center =  ellipse_center(e_obj)

    first_point_angle_to_horizontal = np.arctan2(xyz[0,1] - center[1], xyz[0,0] - center[0]) #angle between the line of center to first point in `xyz` and the horizontal
    # first_angle = np.absolute(first_point_angle_to_horizontal - ellipse_tilt_angle)
    first_angle = first_point_angle_to_horizontal - ellipse_tilt_angle

    second_point_angle_to_horizontal = np.arctan2(xyz[1,1] - center[1], xyz[1,0] - center[0]) #angle between the line of center to second point in `xyz` and the horizontal
    # second_angle = np.absolute(second_point_angle_to_horizontal - ellipse_tilt_angle)
    second_angle = second_point_angle_to_horizontal - ellipse_tilt_angle


    if np.sign(first_angle) == np.sign(second_angle):
        first_angle = np.absolute(first_angle)
        second_angle = np.absolute(second_angle)
        while first_angle > np.pi: first_angle -= np.pi
        while second_angle > np.pi: second_angle -= np.pi

        #At this point both angles are in the range (0, 180)
        # flip if not anticlockwise
        # No, we want to measure angle from the postive x axis
        if first_angle > second_angle: #anti
            first_angle = np.pi - first_angle
        return first_angle
    else:
        first_angle = np.absolute(first_angle)
        if first_angle < np.pi/2: #anti
            first_angle = np.pi - first_angle
        return first_angle




def calculate_eccentricity(traj, smiles = None, calc_convex_hull = True, calc_angle = True):
    indices = _get_largest_ring_indices(traj, smiles)


    # print(indices)
    for frame in traj:
        xyz = frame.xyz[0][indices, :] #only backbone indices

        xyz, variance_ratio = get_pca(xyz)

        e_obj = fit_ellipse(xyz[:,0], xyz[:,1])
        if calc_convex_hull:
            ch_xyz = get_convex_hull(xyz)
            ch_e_obj = fit_ellipse(ch_xyz[:,0], ch_xyz[:,1])

        if calc_angle:
            angle = get_angle_helper(e_obj, xyz)

        if not calc_convex_hull and not calc_angle:
            yield get_eccentricity(e_obj), variance_ratio
        elif not calc_convex_hull and calc_angle:
            yield get_eccentricity(e_obj), variance_ratio, angle
        elif calc_convex_hull and not calc_angle:
            yield get_eccentricity(e_obj), get_eccentricity(ch_e_obj), variance_ratio
        else:
            yield get_eccentricity(e_obj), get_eccentricity(ch_e_obj), variance_ratio, angle

    # except Exception as e: #FIXME ugly
    #     print(e)
    #     pass

def show_ellipse_fitting(md_traj, smiles = None, show_ch = True, save_path = None):
    """
        red colors  the first two points
        the first point has an additional line linked to the center (x) cross
        second cross indicate the center for CH ellipse
        balck dots are convex hull points
        grey are non convex hull point
        dotted ellipse fitted to all points
        solid ellipse fitted to only CH points
    """


    import matplotlib.pyplot as plt
    from ..geometry import plot_ellipse,ellipse_axis_length
    """
        only plots the first frame of md_traj
    """
    traj = md_traj[0]

    indices = _get_largest_ring_indices(traj, smiles)

    xyz = traj.xyz[0][indices, :] #only backbone indices
    # print(xyz)

    xyz, variance_ratio = get_pca(xyz)
    ch_xyz = get_convex_hull(xyz)

    # xyz[:,0], xyz[:,1] = center_2D_points(xyz[:,0], xyz[:,1])
    e_obj = fit_ellipse(xyz[:,0], xyz[:,1])

    # ch_xyz[:,0], ch_xyz[:,1] = center_2D_points(ch_xyz[:,0], ch_xyz[:,1])
    ch_e_obj = fit_ellipse(ch_xyz[:,0], ch_xyz[:,1])



    fig,ax = plt.subplots(1,1, figsize = (10,10))
    plt.hlines(0, -1, 1)
    plt.vlines(0, -1, 1)

    # print(center) #is not exactly at (0,0)
    axes, angle, center = get_info_from_e_obj(e_obj, [ellipse_axis_length, ellipse_angle_of_rotation, ellipse_center])
    plt.scatter(center[0], center[1], c = "black", marker = "x")
    plot_ellipse(semimaj=axes[0],semimin=axes[1],phi=angle, x_cent=center[0], y_cent=center[1], plot_kwargs = {"linestyle" : "-", "color" : "grey"}, ax = ax)
    plot_ellipse(semimaj=axes[0],semimin=axes[1],phi=angle, x_cent=center[0], y_cent=center[1], plot_kwargs = {"linestyle" : "--", "color" : "black"}, ax = ax)

    plt.plot([center[0], xyz[0,0]], [center[1], xyz[0,1]], "grey")

    ################convex hull

    if show_ch:
        axes, angle, center = get_info_from_e_obj(ch_e_obj, [ellipse_axis_length, ellipse_angle_of_rotation, ellipse_center])
        plt.scatter(center[0], center[1], c = "black", marker = "x")
        # plt.scatter(center[0], center[1], c = "grey", marker = "o", fillstyle = "none")

        plot_ellipse(semimaj=axes[0],semimin=axes[1],phi=angle, x_cent=center[0], y_cent=center[1], plot_kwargs = {"linestyle" : "-", "color" : "black"}, ax = ax)
        # plt.plot([center[0], ch_xyz[0,0]], [center[1], ch_xyz[0,1]], "black")


    ################

    plt.scatter(xyz[:,0], xyz[:,1], color = "grey", marker = "o")
    plt.scatter(ch_xyz[:,0], ch_xyz[:,1], color = "black", marker = "o")
    plt.scatter(xyz[0:2,0], xyz[0:2,1], color = "red", marker = "o")


    angle = get_angle_helper(e_obj, xyz)
    if show_ch:
        # Dev_from_3D : {:.2f}
        string = """
        Eccen (roundness): {:.2f}
        Eccen_CH : {:.2f}
        Norm_3rd_PC (flatness): {:.2f}
        Angle : {:.0f}{}
        """.format(get_eccentricity(e_obj), get_eccentricity(ch_e_obj),  3 * variance_ratio[-1], angle*180/3.1415926, u'\N{DEGREE SIGN}')
    else:
        # Dev_from_3D : {:.2f}
        string = """
        Eccen (roundness) : {:.2f}
        Norm_3rd_PC (flatness): {:.2f}
        Angle : {:.0f}{}
        """.format(get_eccentricity(e_obj),  3 * variance_ratio[-1], angle*180/3.1415926, u'\N{DEGREE SIGN}')
    plt.text(-.01,.0, string,
    verticalalignment='bottom',
    horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=20)

    plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False,
    labelleft = False,
    right = False,
    left = False) # labels along the bottom edge are off
    if type(save_path) is str:
        plt.savefig(save_path, transparent = True)
    plt.show()
