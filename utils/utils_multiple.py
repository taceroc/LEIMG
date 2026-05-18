import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import alphashape
from scipy import interpolate
from numba import jit, prange


##utilities.py

def define_bandpass_rubin(wavel):
    """
        Find the bandpass of the wavelenght

        Arguments:
            wavel: wavelenght of LE observation in micrometers
        Returns:
            lsstband: wavelenght limits, band letter, max magnitude, min flux
    """

        
    lsstu = [np.array([3104.9999999999995, 4085.9999999999995]) * 1e-4, 'u', 24.5, 10**((-48.6 - 24.5) / 2.5)]
    lsstg = [np.array([3865.9999999999995, 5669.999999999999]) * 1e-4, 'g', 25, 10**((-48.6 - 25) / 2.5)]
    lsstr = [np.array([5369.999999999999, 7059.999999999999]) * 1e-4, 'r', 25, 10**((-48.6 - 25) / 2.5)]
    lssti = [np.array([6759.999999999999, 8329.999999999998]) * 1e-4, 'i', 24, 10**((-48.6 - 24) / 2.5)]
    lsstz = [np.array([8029.999999999998, 9384.999999999998]) * 1e-4, 'z', 23.5, 10**((-48.6 - 23.5) / 2.5)]
    lssty = [np.array([9083.999999999998, 10944.999999999998]) * 1e-4, 'y', 23.5, 10**((-48.6 - 23.5) / 2.5)]

    for bandpasses in [lsstu, lsstg, lsstr, lssti, lsstz, lssty]:
        if bandpasses[0][0] <= wavel <= bandpasses[0][1]:
            return bandpasses


def project_skyplane_plane_cube_2ddust(x_inter, y_inter, z_inter, cossigma):
    """
    Project x,y,z of intersection paraboloid and dust into the sky plane
    
    Arguments:
        x_inter, y_inter, z_inter: intersection paraboloid + dust in ly
        cossigma: cos scattering angle
    
    Return:
        x_img_t, y_img_t, z_img_t:  position sky plane
    """
    norm_r = np.sqrt(x_inter**2 + y_inter**2 + z_inter**2)
    xdotr = x_inter / norm_r #phi angle
    phis = np.arccos(xdotr)#, cossigma >> scatt >> 
    scatt = np.arccos(cossigma)
    # transform to sky plane, RR patat paper
    x_img_t = (x_inter * np.cos(scatt) * np.cos(phis) 
               + y_inter * np.sin(scatt) * np.cos(phis)
               + z_inter * np.sin(phis))
    y_img_t = -x_inter * np.sin(phis) + y_inter * np.cos(phis)
    z_img_t = (- x_inter * np.sin(scatt) * np.cos(phis) 
               - y_inter * np.sin(scatt) * np.sin(phis) 
               + z_inter * np.cos(scatt)) 
    
    return x_img_t, y_img_t, z_img_t

def convert_ly_to_arcsec(d, r_le_in):
    import astropy.cosmology.units as cu
    from astropy import units as units
    units.add_enabled_units(cu)
    
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
    d = (d * units.lyr).to(units.Mpc)
    reds = d.to(cu.redshift, cu.redshift_distance(cosmo, kind="comoving"))
    # linear size = angular_size * d_A
    d_A = cosmo.angular_diameter_distance(reds)
    r_le_in = (r_le_in * units.lyr).to(units.Mpc)
    r_le_in = (r_le_in / d_A ).value / (np.pi / 180 / 3600)
    return r_le_in

def convert_arcsec_to_ly(d, r_le_in):
    import astropy.cosmology.units as cu
    from astropy import units as units
    units.add_enabled_units(cu)
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
    d = (d * units.lyr).to(units.Mpc)
    reds = d.to(cu.redshift, cu.redshift_distance(cosmo, kind="comoving"))
    # linear size = angular_size * d_A
    d_A = cosmo.angular_diameter_distance(reds)
    r_le_in = r_le_in * (np.pi / 180 / 3600) * d_A.value 
    r_le_in = (r_le_in * units.Mpc).to(units.lyr)
    return r_le_in

def define_shape_to_interpolate(x_inter_arcsec_total, y_inter_arcsec_total, pixel = 0.2):

    x_lim_min, x_lim_max = np.min(x_inter_arcsec_total) , np.max(x_inter_arcsec_total)
    y_lim_min, y_lim_max = np.min(y_inter_arcsec_total), np.max(y_inter_arcsec_total)
    
    x_tot_arcsec = np.abs(int(x_lim_max - x_lim_min))
    y_tot_arcsec = np.abs(int(y_lim_max - y_lim_min))
    
    x_size_img = int(x_tot_arcsec / pixel)
    y_size_img = int(y_tot_arcsec / pixel)
    print(x_size_img, y_size_img)

    x_all, y_all = np.meshgrid(np.linspace(x_lim_min, x_lim_max, y_size_img),
                           np.linspace(y_lim_min, y_lim_max, x_size_img ))

    points = list(zip(x_inter_arcsec_total, y_inter_arcsec_total))
    index = [i[0] for i in sorted(enumerate(points), key=lambda x: [x[1][1], x[1][0]])]
    arr_points = np.array(points)
    mpoints = [(X, Y) for X, Y in list(zip(arr_points[index, 0], arr_points[index, 1]))]
    alpha_shape = alphashape.alphashape(mpoints, alpha=0.2)

    return alpha_shape, x_size_img, y_size_img, arr_points, index, x_all, y_all

def interpolation_surface(new_xs, new_ys, surface_original):
    surface_inter_y = interpolate.LinearNDInterpolator(list(zip(new_xs, new_ys)), surface_original)
    # surface_inter_y = interpolate.NearestNDInterpolator(list(zip(new_xs, new_ys)), surface_original)
    return surface_inter_y

@jit(nopython=True)
def cal_inter(x_p, y_p, z_p, h, k, l, ct, r0ly, dt0):
    intersection_points = []
    for i in range(len(x_p)):
        for j in range(len(y_p)):
            x_par, y_par, z_par = x_p[i,j], y_p[i,j], z_p[i,j]

            # Check if the point is inside both the sphere and the paraboloid
            sphere_condition = (r0ly - 0.05)**2 <= ((x_par - h)**2 + (y_par - k)**2 + (z_par - l)**2) <= (r0ly + 0.05)**2
            paraboloid_condition = ((ct - dt0 )**2 + 2 * (ct - dt0) * z_par) <= (x_par**2 + y_par**2) <= (ct**2 + 2 * ct * z_par)

            if (sphere_condition and paraboloid_condition):
                intersection_points.append((x_par, y_par, z_par))
    return intersection_points


@jit(nopython=True)
def cal_inter_cube(x_p, y_p, z_p, x_min, x_max, y_min, y_max, z_min, z_max, ct, dt0):
    intersection_points = []
    for i in range(len(x_p)):
        for j in range(len(y_p)):
            x_par, y_par, z_par = x_p[i,j], y_p[i,j], z_p[i,j]

            # Check if the point is inside both the sphere and the paraboloid
            cube_condition = (x_min<=x_par<=x_max) and (y_min<=y_par<=y_max) and (z_min<=z_par<=z_max)
            paraboloid_condition =  ((ct - dt0 )**2 + 2 * (ct - dt0) * z_par) <= (x_par**2 + y_par**2) <= (ct**2 + 2 * ct * z_par)

            if (cube_condition and paraboloid_condition):
                intersection_points.append((x_par, y_par, z_par))
    return intersection_points


@jit(nopython=True)
def cal_inter_sheetdust(x_all, y_all, z_all, r_in, r_out, act, bct, x_min, x_max, y_min, y_max, z_min, z_max, params):
    x_all_inside = []
    y_all_inside = []
    z_all_inside = []
    # print(r_in.shape)
    for i in range(x_all.shape[-1]):
        for j in range(y_all.shape[-1]):
            if (r_in[j,i] <= np.sqrt((x_all[j,i] + act)**2 + (y_all[j,i] + bct)**2) <= r_out[j,i]):
                if (-1E-5 <= z_all[j,i] * params[2] + params[-1] + params[0]*x_all[j,i] + params[1]*y_all[j,i] <= 1E-5):
                    if (x_min<=x_all[j,i]<=x_max) and (y_min<=y_all[j,i]<=y_max) and (z_min<=z_all[j,i]<=z_max):
                        x_all_inside.append(x_all[j,i])
                        y_all_inside.append(y_all[j,i])
                        z_all_inside.append(z_all[j,i])
    return x_all_inside, y_all_inside, z_all_inside



def bin_data_xyz(x,y,z):
    #find the unique x,y coordiantes
    pairs = np.vstack([x, y]).T
    unique_pairs, ind, inverse, count = np.unique(pairs, axis=0, return_counts=True, return_index=True, return_inverse=True)
    # count = 44*int((count - np.min(count)) / (np.max(count)-np.min(count))
    # count = 44*np.array((count - count.min()) / (count.max() - count.min())).astype(int) +1
    print(unique_pairs.shape, count.max())
    if count.max() > 5:
        size_depth = count.max()
    else:
        size_depth = 10
    #create an initial array the same size as the number of unique x,y values and the number of times that each pair appear
    arr_init = np.zeros([len(unique_pairs),size_depth])
    ratio = int((abs(y.min() - y.max()) / abs(x.min() - x.max())))
    print("ratio", abs(x.min() - x.max()))
    size_x = int(np.sqrt(arr_init.shape[0]/2))
    if ratio != 0:
        size_y = ratio*size_x
    else:
        size_y = size_x
    print((size_x * size_y) <= len(unique_pairs))
    #index the x,y,z values
    #create arrays with a size of int(sqrt(len unique arrays))
    x_bins = np.linspace(x.min(), x.max(), size_x)
    y_bins = np.linspace(y.min(), y.max(), size_y)

    #what a great function, return an array of size=original and the values are the bin at which that value belong
    x_indices = np.digitize(x, x_bins)
    y_indices = np.digitize(y, y_bins)

    #same for z, but of size count
    z_bins = np.linspace(z.min(), z.max(), size_depth)
    z_indices = np.digitize(z, z_bins)
    print(np.unique(z_indices))

    #ind: size of the unique arrays: the values are the index on the original array
    #inverse: size of the original array: the values are the index on the unique array
    #arr_init: contains the z values for each x,y unique pair. the depth is count and it is organized, each of the 11 is a unique value of z
    for i, val in enumerate(ind):
        if count[i] > 1:
            for j in range(count[i]-1):
                arr_init[i,z_indices[np.where(inverse == i)[0][j]]-1] = z_bins[z_indices[np.where(inverse == i)[0][j]]-1]
        else:
            # print(i,val,z[val])
            arr_init[i,z_indices[val]-1] = z_bins[z_indices[val]-1]
    # replace 0 for nan to extract later the min !=0
    arr_init[arr_init == 0] = 'nan'

    return x_indices, y_indices, arr_init, unique_pairs, inverse, count, size_x, size_y


def find_mask_angles(ini_angle, end_angle, surface_val, x_img, y_img):

    def find_quadrant(angle):
        div = angle//90
        if div == 0:
            return 1
        elif div == 1:
            return 2
        elif div == 2:
            return 3
        elif div >= 3:
            return 4

   

    def find_limits(ini_angle, end_angle, x_img, y_img):
        xs_ini = np.abs(x_img.min() - x_img.max())*np.cos(np.deg2rad(ini_angle))/2 #rs_pc*np.cos(np.deg2rad(ini_angle)) + ((abs(le_img.new_xs.max()) - abs(le_img.new_xs.min()))/2)
        ys_ini = np.abs(y_img.min() - y_img.max())*np.sin(np.deg2rad(ini_angle))/2 #rs_pc*np.sin(np.deg2rad(ini_angle)) + ((abs(le_img.new_ys.max()) - abs(le_img.new_ys.min()))/2)
        xs_end = np.abs(x_img.min() - x_img.max())*np.cos(np.deg2rad(end_angle))/2 #rs_pc*np.cos(np.deg2rad(end_angle)) + ((abs(le_img.new_xs.max()) - abs(le_img.new_xs.min()))/2)
        ys_end = np.abs(y_img.min() - y_img.max())*np.sin(np.deg2rad(end_angle))/2 #rs_pc*np.sin(np.deg2rad(end_angle)) + ((abs(le_img.new_ys.max()) - abs(le_img.new_ys.min()))/2)

        xs_min = (x_img.min() + np.abs(x_img.min() - x_img.max())/2) + np.min([xs_ini, xs_end])
        xs_max = (x_img.max() - np.abs(x_img.min() - x_img.max())/2) + np.max([xs_ini, xs_end])

        if xs_min <= x_img.min()+np.abs(x_img.min() - x_img.max())/2:
            if xs_max <= x_img.min()+np.abs(x_img.min() - x_img.max())/2:
                xs_max = x_img.min()+np.abs(x_img.min() - x_img.max())/2
        if xs_max >= x_img.min()+np.abs(x_img.min() - x_img.max())/2:
            if xs_min >= x_img.min()+np.abs(x_img.min() - x_img.max())/2:
                xs_min = x_img.min()+np.abs(x_img.min() - x_img.max())/2
        
        ys_min = (y_img.min() + np.abs(y_img.min() - y_img.max())/2) + np.min([ys_ini, ys_end])
        ys_max = (y_img.max() - np.abs(x_img.min() - x_img.max())/2) + np.max([ys_ini, ys_end])
        
        
        x_mask = ((x_img >= xs_min)&(x_img <= xs_max))
        y_mask = ((y_img >= ys_min)&(y_img  <= ys_max))
        
        # surface_val = surface_val * x_mask*y_mask

        return x_mask*y_mask

    def find_limits_quadrant(x_img, y_img, quadrant, angle=[], last=False, first=False):
        if ((last == False) & (first == False)):
            print(f'MID QUADRANT')
            if quadrant == 1:
                ini, end = 0, 89.9
            if quadrant == 2:
                ini, end = 90, 179.9
            if quadrant == 3:
                ini, end = 180, 269.9
            if quadrant == 4:
                ini, end = 270, 359.9
        elif last==True:
            print(f'LAST QUADRANT')
            if quadrant == 1:
                ini, end = 0, angle[1]
            if quadrant == 2:
                ini, end = 90, angle[1]
            if quadrant == 3:
                ini, end = 180, angle[1]
            if quadrant == 4:
                ini, end = 270, angle[1]
        elif first==True:
            print(f'FIRST QUADRANT')
            if quadrant == 1:
                ini, end = angle[0], 89.9
            if quadrant == 2:
                ini, end = angle[0], 179.9
            if quadrant == 3:
                ini, end = angle[0], 269.9
            if quadrant == 4:
                ini, end = angle[0], 359.9
        return find_limits(ini, end, x_img, y_img)

    
    if (ini_angle == 0) & (end_angle == 360):
        xs_min = x_img.min()
        xs_max = x_img.max()

        ys_min = y_img.min()
        ys_max = y_img.max()

        x_mask = ((x_img >= xs_min)&(x_img <= xs_max))
        y_mask = ((y_img >= ys_min)&(y_img  <= ys_max))

        mask = x_mask*y_mask
            
    else:
        if end_angle < 0:
            end_angle = end_angle+360
        mask = np.zeros_like(surface_val)
        ini_quadrant = find_quadrant(ini_angle)
        end_quadrant = find_quadrant(end_angle)
        total_quadrants = np.arange(ini_quadrant, end_quadrant+1, 1)
        if len(total_quadrants) == 1:
            mask = find_limits(ini_angle, end_angle, x_img, y_img)
        else:
            for i in total_quadrants:
                print(f"QUADRANT {i}: {ini_angle}, {end_angle}")
                if i == total_quadrants[0]:
                    print(f'FOR FIRST QUADRANT')
                    mask = mask + find_limits_quadrant(x_img, y_img, i, angle=[ini_angle,end_angle], last=False, first=True)
                elif i == total_quadrants[-1]:
                    print(f'FOR LAST QUADRANT')
                    mask = mask + find_limits_quadrant(x_img, y_img, i, angle=[ini_angle,end_angle], last=True, first=False)
                else:
                    print(f'FOR MID QUADRANT')
                    mask = mask + find_limits_quadrant(x_img, y_img, i, angle=[ini_angle,end_angle], last=False, first=False)

    mask[mask >= 1] = 1
    if end_angle < 0:
        mask = ~mask

    return mask
    