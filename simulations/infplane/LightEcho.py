import sys
import logging
import numpy as np
from utils import fix_constants as fc

logger = logging.getLogger(__name__)

class LE:
    """
        LE define the properties of the Light Echo
        Subclasses:
            LEPlane
            LESphereCentered
            LESphericalBulb*** incomplete
            LESheetDust
    """
    def __init__(self, ct, source):
        """
            x: initialize x positions in ly
            ct: time of LE observation in years
        """
        self.x = 0
        self.ct = ct
        self.r_le2 = 0
        self.x_inter_values = 0
        self.y_inter_values = 0
        self.z_inter_values = 0
        self.x_projected = 0
        self.y_projected = 0
        self.dt0 = source.dt0
        self.d = source.d   
        self.theta = [0, 2*np.pi]
        
    def calculate_intersection_x_y(self, y_1):
        """
            Calculate the intersection points x,y between a DustShape and the paraboloid
        """
        y_2 = -1*y_1
        # keep no nan values
        y_inter = np.hstack((y_1, y_2))
        self.y_inter_values = y_inter[~np.isnan(y_inter)]
        # extract x where y is no nan
        x_inv_nan = np.hstack((self.x, self.x.copy()))
        self.x_inter_values = x_inv_nan[~np.isnan(y_inter)]

        
    def calculate_rho_thickness(self):
        """
            Calculate the rho in the sky plane (xy) and thickness according to Sugermann 2003 and Xi 1994
            and calculate the radius out and radius min according to that thickness
        """
        r_le = np.sqrt(self.r_le2)
        rhos = np.sqrt(2 * self.z_inter_values * self.ct + (self.ct)**2 )
        half_obs_thickness = np.sqrt( (self.ct / rhos) ** 2 * self.dz0 ** 2 ) / 2 #+ ( (rhos * fc.c / (2 * self.ct)) + ( fc.c * self.ct / (2 * rhos) )) ** 2 * self.dt0  ** 2 ) / 2
        
        # -- include the thickness in xy plane
        if len(np.array(half_obs_thickness)) != 1:
            self.r_le_out = np.ones(len(half_obs_thickness)) * r_le 
            self.r_le_in = r_le - 2*half_obs_thickness
        else:
            self.r_le_out = r_le 
            self.r_le_in = r_le - 2*half_obs_thickness
        
        
    def final_xy_projected(self):
        """
        Calculate the x,y points in arcseconds
        Only valid when the dust and the paraboloid have a analytical expresion (and the analtyical expression is a circumference)

        Arguments:
            phis: angle in the sky plane
            r_le_out, r_le_in: out and inner radii in arcsec
            act: center of LE in arcsec

        Returns:
            x_projected, y_projected: x,y position in the x-y plane in arcseconds
        """
        self.calculate_rho_thickness()
        phis = np.arctan2(self.y_inter_values, self.x_inter_values)     
        radii_p = [self.r_le_out, self.r_le_in]
        # self.r_le = np.ones_like(self.y_inter_values)*radii_p
        if self.F != 0:
            coefx = (self.A/self.F)*self.ct
            coefy = (self.B/self.F)*self.ct
        else:
            coefx = 0
            coefy = 0
        xs_p = np.concatenate([radii_p[0] * np.cos(phis) - coefx, radii_p[1] * np.cos(phis) - coefx]).reshape(2, len(phis))
        ys_p = np.concatenate([radii_p[0] * np.sin(phis) - coefy, radii_p[1] * np.sin(phis) - coefy]).reshape(2, len(phis))

        self.x_projected = xs_p.reshape(1,2,len(phis))
        self.y_projected = ys_p.reshape(1,2,len(phis))


        return self.x_projected, self.y_projected
    
    def run(self):
        """
            Return x,y,z intersection in ly and x,y projected in the sky plane in arcseconds
        """
        self.x_inter_values, self.y_inter_values, self.z_inter_values = self.get_intersection_xyz()
        # logger.info(f"There are %s intersections points in x", self.x_inter_values.shape)
        # logger.info(f"There are %s intersections points in y", self.y_inter_values.shape)
        # logger.info(f"There are %s intersections points in z", self.z_inter_values.shape)

        def process_list(my_list):
            if my_list[0] == 0:
                raise ValueError("No intersection points. Cannot continue.")
            print("Intersections")
        try:
            process_list(self.x_inter_values.shape)
        except ValueError as e:
            print(e)
            return e
        self.x_projected, self.y_projected = self.final_xy_projected()

        self.z_projected = self.func_for_z(self.x_projected, self.y_projected) #-(self.D/self.F) - (self.A/self.F) * self.x_inter_values - (self.B/self.F) * self.y_inter_values
        
        return self.x_inter_values, self.y_inter_values, self.z_inter_values, self.x_projected, self.y_projected, self.z_projected
    

class LEPlane(LE):
    """
        LE_plane(LE) defines a subclass of LE
            LE_plane
    """
    def __init__(self, ct, plane, source):
        """
            x: initialize x positions in ly
            ct: time of LE observation in years
            plane.eq_params = [A, B, C, D], Ax + By + Fz + D = 0
            plane.dz0: depth inf plane of dust in ly
            
        """
        super().__init__(ct, source)
        self.A = plane.eq_params[0]
        self.B = plane.eq_params[1]
        self.F = plane.eq_params[2]
        self.D = plane.eq_params[3]
        self.dz0 = plane.dz0
        self.r_le_out = 0
        self.r_le_in = 0
        self.check_time_forLE()
        self.func_for_z = self.func_for_z_plane
        self.x_min_lim = None
        self.x_max_lim = None
        self.y_min_lim = None
        self.y_max_lim = None
        self.z_min_lim = None
        self.z_max_lim = None

    def check_time_forLE(self):
        """
            If the plane is behind the source check that the ct time is later than the start of the LE
        """
        if not isinstance(self.D, (np.ndarray, list)):
            logger.info(True)
            if -(self.D/self.F) < 0:
                ti = (2 * (self.D/self.F))/(fc.c * (1 + (self.A/self.F)**2 + (self.B/self.F)**2))
                if ti >= self.ct:
                    logger.warning(f"There is no LE at {self.ct} years, LE starts to expand at {ti/fc.dtoy} days or {ti} years")
                    raise ValueError('No LE yet')
        
    def calculate_rle2(self):
        """
            Calcualte the radii square of the resultant LE plane
        """
        self.r_le2 = -2 * self.ct * (self.D/self.F) + (self.ct)**2 * (1 + (self.B/self.F)**2 + (self.A/self.F)**2)
        return self.r_le2
    
    def get_intersection_xyz(self):
        """
            Calculate the intersection points x,y,z between a DustShape and the paraboloid
        """
        self.calculate_rle2()
        theta_p = np.linspace(self.theta[0], self.theta[1], 1000)
        logger.info('Angle to create intersection')
        logger.info(f'{self.theta[0]}, {self.theta[1]}')
     
        
        self.x_inter_values = np.sqrt(self.r_le2) * np.cos(theta_p) - (self.A/self.F)*self.ct
        self.y_inter_values = np.sqrt(self.r_le2) * np.sin(theta_p) - (self.B/self.F)*self.ct
        # calculate z = z0 - ax >> plane equation
        self.z_inter_values = self.func_for_z_plane(self.x_inter_values, self.y_inter_values) #-(self.D/self.F) - (self.A/self.F) * self.x_inter_values - (self.B/self.F) * self.y_inter_values
        
        return self.x_inter_values, self.y_inter_values, self.z_inter_values
    
    def func_for_z_plane(self, val_x=0, val_y=0):
        return -(self.D/self.F) - (self.A/self.F) * val_x - (self.B/self.F) * val_y
      