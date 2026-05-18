import numpy as np
import logging
from scipy.interpolate import griddata
import utils.utils_multiple as utils_multiple

logger = logging.getLogger(__name__)

class SurfaceAssembler:
    def __init__(self, x_bins: np.ndarray, y_bins: np.ndarray, all_surface: np.ndarray, surface_known: list):
        self.x_bins = x_bins
        self.y_bins = y_bins
        self.all_surface = all_surface
        self.surface_known = surface_known

    def create_le_surface_interpolate(self, x_known, y_known): #in le_img outside loop
        """
            Interpolate the surface values to fill the image, return also the x,y  values in arc
            Arguments:
              
            Return:
                inter_surface: the light echo surface image in units ly, cm,
                x_img: x axis in arc
                y_img: y axis in arc
                z_img_ly: z in ly
                
        """
        points = np.vstack([np.array(x_known), np.array(y_known)]).T
        
        grid_x, grid_y = np.mgrid[
        self.x_bins.min():self.x_bins.max():self.x_bins.shape[0]*1j, 
        self.y_bins.min():self.y_bins.max():self.y_bins.shape[0]*1j
        ]
        logger.info('griddata %s',(np.array(self.all_surface).shape, points.shape, grid_x.shape, grid_y.shape ))

        # surface_known = np.concatenate([self.all_surface, self.all_surface])
        if len(self.surface_known) == points.shape[0]:
            values_for_grid = np.array(self.surface_known)
            logger.info('know')
            logger.info(f"{values_for_grid.min()}, {values_for_grid.max()}")
        else:
            logger.info('all')
            values_for_grid = np.array(self.all_surface)
            logger.info(f"{values_for_grid.min()}, {values_for_grid.max()}")
            

        logger.info(f'griddata %s', (values_for_grid.shape, points.shape, grid_x.shape, grid_y.shape) )
        logger.info(f"{values_for_grid.min()}, {values_for_grid.max()}")

        inter_surface = griddata(points, values_for_grid, (grid_x.T, grid_y.T), method='nearest')
        
        print('before masking too big')
        print(inter_surface.min(), inter_surface.max())

        return inter_surface, grid_x, grid_y

    def apply_annulus_mask(self, inter_surface: np.ndarray, grid_x: np.ndarray, grid_y: np.ndarray, 
                           r_in: float, r_out: float, center_x: float, center_y: float, z_all_ly) -> np.ndarray:

        # center_x = np.nanmean(abs(self.results_simulation.act_all[0]-self.results_simulation.act_all[-1]))
        # center_y = np.nanmean(abs(self.results_simulation.bct_all[0]-self.results_simulation.bct_all[-1]))
        
        distance_from_center = np.sqrt((grid_x.T - center_x)**2 + (grid_y.T - center_y)**2)
        
        # Robust radii from valid values only
        # r_in = np.nanmedian(np.asarray(self.results_simulation.all_r_le_in_arc[-1], dtype=float))
        # r_out = np.nanmedian(np.asarray(self.results_simulation.all_r_le_out_arc[1], dtype=float))

        annulus_mask = (distance_from_center >= r_in) & (distance_from_center <= r_out)


        inter_surface[~annulus_mask] = np.nan

        logger.info('center x,y: %s', (center_x, center_y))
        logger.info("r_in, r_out: %s", (r_in, r_out))
        logger.info("distance min/max: %s", (np.nanmin(distance_from_center), np.nanmax(distance_from_center)))
        logger.info("annulus kept fraction: %s", np.mean(annulus_mask))

        mask_binary = np.isnan(inter_surface)

        logger.info('after masking too big')
        logger.info(f"{np.nanmin(inter_surface)}, {np.nanmax(inter_surface)}")

        z_grid_ly = z_all_ly
        if z_grid_ly.shape != mask_binary.shape:
            if z_grid_ly.T.shape == mask_binary.shape:
                z_grid_ly = z_grid_ly.T
            else:
                raise ValueError(
                    f"z-grid shape {z_grid_ly.shape} does not match surface mask shape {mask_binary.shape}"
                )

        return inter_surface, grid_x.T*(~mask_binary), grid_y.T*(~mask_binary), z_grid_ly*(~mask_binary)
    
    def apply_angle_mask(self, surface_val: np.ndarray, x_img_arc: np.ndarray, y_img_arc: np.ndarray, ini_angle: float, end_angle: float) -> np.ndarray:
        """
            Crop the LE image to the given angular size (the arc)
            Arguments:
                surface_val: the light echo surface image in units ly, cm,
                x_img: x axis in arc
                y_img: y axis in arc
                
            Return:
                surface_val: with the angular size
                
        """

        surface_val_copy = surface_val.copy()

        mask = utils_multiple.find_mask_angles(ini_angle, end_angle, surface_val, x_img_arc, y_img_arc)
        logger.info('ANGLES')
        logger.info(f"{ini_angle}, {end_angle}")

        # ONLY END ANGLE CAN BE NEGATIVE
        if end_angle < 0:
            surface_val = surface_val*(~np.array(mask).astype(bool))
        else:
            surface_val = surface_val*mask
        logger.info('shape surface image final')
        logger.info(f"{surface_val.shape}")
    
        if np.all(mask == False):
            surface_val = surface_val_copy
    
        return surface_val