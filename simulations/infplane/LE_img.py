import sys
import numpy as np

import utils.utils_multiple as utils_multiple

from .DustShape import SphericalBlub, PlaneDust
import logging
logger = logging.getLogger(__name__)

class LEImage:
    """ 
        Initialize object to create LE image
    """
    def __init__(self, LE_geometryanalyticalsource, surface, pixel_resolution=0.2):
        self.pixel = pixel_resolution # arcsec
        self.LE_geom = LE_geometryanalyticalsource
        self.new_xs = utils_multiple.convert_ly_to_arcsec((LE_geometryanalyticalsource.d)+LE_geometryanalyticalsource.z_projected, LE_geometryanalyticalsource.x_projected)
        self.new_ys = utils_multiple.convert_ly_to_arcsec((LE_geometryanalyticalsource.d)+LE_geometryanalyticalsource.z_projected, LE_geometryanalyticalsource.y_projected)

        self.surface_original = surface
        self.surface_val = 0
        self.surface_img = 0
        
    def _estimate_image_size(self):

        x_lim_min, x_lim_max = np.min(self.new_xs), np.max(self.new_xs)
        y_lim_min, y_lim_max = np.min(self.new_ys), np.max(self.new_ys)

        x_tot_arcsec = np.abs(round((x_lim_max - x_lim_min),0))
        y_tot_arcsec = np.abs(round((y_lim_max - y_lim_min),0))
        logger.info("size arcsec %s", (x_tot_arcsec, y_tot_arcsec))
        
        x_size_img = int(x_tot_arcsec / self.pixel)
        y_size_img = int(y_tot_arcsec / self.pixel)
        logger.info("size img pixels %s", (x_size_img, y_size_img))

        return x_size_img, y_size_img
    
    def create_image_grid(self, x_size_img, y_size_img):
        x_lim_min, x_lim_max = np.min(self.new_xs), np.max(self.new_xs)
        y_lim_min, y_lim_max = np.min(self.new_ys), np.max(self.new_ys)

        x_lim_min_ly, x_lim_max_ly = np.min(self.LE_geom.x_projected), np.max(self.LE_geom.x_projected)
        y_lim_min_ly, y_lim_max_ly = np.min(self.LE_geom.y_projected), np.max(self.LE_geom.y_projected)
        
        
        x_all, y_all = np.meshgrid(np.linspace(x_lim_min, x_lim_max, x_size_img),
                                   np.linspace(y_lim_min, y_lim_max, y_size_img ))
        
        x_bins = np.linspace(x_lim_min, x_lim_max, x_size_img)
        y_bins = np.linspace(y_lim_min, y_lim_max, y_size_img)
        logger.info('min max for phase 0 x%s', (x_lim_min, x_lim_max))
        logger.info('min max for phase 0 y%s', (y_lim_min, y_lim_max))
        

        x_all_ly, y_all_ly = np.meshgrid(np.linspace(x_lim_min_ly, x_lim_max_ly, x_size_img), 
                                         np.linspace(y_lim_min_ly, y_lim_max_ly, y_size_img ))
        z_all_ly = self.LE_geom.func_for_z(x_all_ly, y_all_ly)

        x_max, x_min = np.max(self.LE_geom.x_inter_values[(self.new_xs[0,0,:] <= x_lim_max) & 
            (self.new_xs[0,1,:] >= x_lim_min)]), np.min(self.LE_geom.x_inter_values[(self.new_xs[0,0,:] <= x_lim_max) & 
            (self.new_xs[0,1,:] >= x_lim_min)])
        y_max, y_min = np.max(self.LE_geom.y_inter_values[(self.new_ys[0,0,:] <= y_lim_max) & (self.new_ys[0,1,:] >= y_lim_min)]), np.min(self.LE_geom.y_inter_values[(self.new_ys[0,0,:] <= y_lim_max) & (self.new_ys[0,1,:] >= y_lim_min)])

        x_max_p, x_min_p = np.max(self.LE_geom.x_projected[(self.new_xs <= x_lim_max) & (self.new_xs >= x_lim_min)]), np.min(self.LE_geom.x_projected[(self.new_xs <= x_lim_max) & (self.new_xs >= x_lim_min)])
        y_max_p, y_min_p = np.max(self.LE_geom.y_projected[(self.new_ys <= y_lim_max) & (self.new_ys >= y_lim_min)]), np.min(self.LE_geom.y_projected[(self.new_ys <= y_lim_max) & (self.new_ys >= y_lim_min)])

        self.new_xs = self.new_xs[(self.new_xs <= x_lim_max) & (self.new_xs >= x_lim_min)]
        self.new_ys = self.new_ys[(self.new_ys <= y_lim_max) & (self.new_ys >= y_lim_min)]

        return x_bins, y_bins, x_all.T, y_all.T, x_all_ly.T, y_all_ly.T, z_all_ly.T, [x_max, x_min], [y_max, y_min], [x_max_p, x_min_p], [y_max_p, y_min_p]


class LEImageAnalytical(LEImage):
    """
        Subclass when solution is analytical: infinite plane and sphere centered
    """

    def __init__(self, LE_geometryanalyticalsource, geometry, surface, pixel_resolution = 0.2):
        super().__init__(LE_geometryanalyticalsource, surface, pixel_resolution)
        # inner an outer radii of LE
        self.r_le_in = utils_multiple.convert_ly_to_arcsec(LE_geometryanalyticalsource.d, LE_geometryanalyticalsource.r_le_in)
        self.r_le_out = utils_multiple.convert_ly_to_arcsec(LE_geometryanalyticalsource.d, LE_geometryanalyticalsource.r_le_out)

        self.geometry_to_use = geometry
        # act, bct: origin of LE in x and y
        self.act = LE_geometryanalyticalsource.ct * (self.geometry_to_use.eq_params[0] / self.geometry_to_use.eq_params[2]) if self.geometry_to_use.eq_params[2] != 0 else 0
        self.act = utils_multiple.convert_ly_to_arcsec(LE_geometryanalyticalsource.d, self.act)
        self.bct = LE_geometryanalyticalsource.ct * (self.geometry_to_use.eq_params[1] / self.geometry_to_use.eq_params[2]) if self.geometry_to_use.eq_params[2] != 0 else 0
        self.bct = utils_multiple.convert_ly_to_arcsec(LE_geometryanalyticalsource.d, self.bct)

    
