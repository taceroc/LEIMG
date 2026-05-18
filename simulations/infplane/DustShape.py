import numpy as np

class DustShape:
    """
        DustShape define the geometrical characteristics of the dust
        Subclasses:
            Infinite Plane
            Sphere centered on the source
    """
    
    def __init__(self, eq_params, dz0):
        # self.check_plane(eq_params)
        self.check_size(eq_params, 4)
        self.eq_params = eq_params
        self.dz0 = dz0
        
    def check_size(self, var, size = 4):
        """
            Check that the arrays have the correct size and type
        """
        if not isinstance(var, (np.ndarray, list)):
            print('Enter a list or an array')
            raise ValueError('Enter a list or an array')
        elif len(var) != size:
            print(f'Array or list {var} should contain {size} parameters')
            raise ValueError('Incorrect size')

class InfPlane(DustShape):
    """
        InfPlane(DustShape) define a subclass of DustShape
            Infinite Plane
    """
    def __init__(self, eq_params, dz0):
        """
            eq_params = [A, B, C, D], Ax + By + Fz + D = 0
            dz0: depth inf plane of dust in ly
        """
        super().__init__(eq_params, dz0)
        self.check_plane(self.eq_params)
        
    def check_plane(self, var):
        """
            Check that the plane is not parallel to z
        """
        if ((int(var[2]) == 0)):
            print('Plane cannot be parallel to z')
            raise ValueError('Plane cannot be parallel to z')
                        
class SphereCenter(DustShape):
    """
        SphereCenter(DustShape) define a subclass of DustShape
            Sphere centered on the source
    """
    def __init__(self, eq_params, dz0):
        """
            eq_params = [A, B, C, D], (x-A)2 + (y-B)2 + (z-F)2 = D^2
            dz0: depth sphere of dust in ly
        """
        super().__init__(eq_params, dz0)
        self.check_centered(self.eq_params)
    
    def check_centered(self, var):
        """
            Check that the sphere is centered at the source
        """
        if ((int(var[0]) != 0) or (int(var[1]) != 0) or (int(var[2]) != 0)):
            print('Sphere of dust most be centered at the source. Try class BLA if the sphere is not centered at the source')
            raise ValueError('Not centered')

class SphericalBlub(DustShape):
    
    def __init__(self, eq_params, dz0):
        """
            eq_params = [A, B, C, D], (x-A)2 + (y-B)2 + (z-F)2 = D^2
            dz0: depth sphere of dust in ly
        """
        super().__init__(eq_params, dz0)


class PlaneDust(DustShape):
    
    def __init__(self, eq_params, dz0, dust_shape, dust_position, size_dustsheet):
        """
            Defines/calculate the limits in x and y direction of the sheet of dust
            
            eq_params = [A, B, C, D], Ax + By + Fz + D = 0
            dz0: depth sphere of dust in ly
            dust_shape: pixel size of the dust (image extracted from Spitzer)
            dust_position: xyz position of the center of the dust cube from the source
            size_dustsheet: size xy of the sheet in ly
        """
        super().__init__(eq_params, dz0)
        super().check_size(dust_position, size = 3)
        super().check_size(size_dustsheet, size = 2)
        self.side_x = dust_shape[0]
        self.side_y = dust_shape[1]
        self.zdepths = dz0
        # if Ax = By = 0, plane perpendicular to z 
        if eq_params[0] == eq_params[1] == 0:
            self.x1 = dust_position[0] + size_dustsheet[0]/2
            self.x2 = dust_position[0] - size_dustsheet[0]/2
            self.y1 = dust_position[1] + size_dustsheet[1]/2
            self.y2 = dust_position[1] - size_dustsheet[1]/2
            self.z12 = dust_position[2]
            self.x_min = min(self.x1, self.x2)
            self.x_max = max(self.x1, self.x2)
            self.y_min = min(self.y1, self.y2)
            self.y_max = max(self.y1, self.y2)
            self.z_min = self.z12
            self.z_max = self.z12
        # if plane is not perpendicular, calculate 4 points that form the dust sheet
        else:
            # assume v1=(A, -B, 0)
            u1n = np.sqrt(eq_params[0]**2 + eq_params[1]**2)
            # normalize v1
            u1 = np.array([eq_params[1] / u1n, -eq_params[0] / u1n, 0 / u1n])
            # another vector perpendicular to v1 and the normal, v2 = nxv1
            u2n = np.sqrt((eq_params[0]*eq_params[2])**2 + (eq_params[1]*eq_params[2])**2 + (-eq_params[0]**2 - eq_params[1]**2)**2)
            u2 = np.array([(eq_params[0]*eq_params[2]) / u2n, (eq_params[1]*eq_params[2]) / u2n, (-eq_params[0]**2 - eq_params[1]**2) / u2n])
            # points are at a distance=size dust sheet over sqrt(2)
            p1 = dust_position + size_dustsheet[0]/np.sqrt(2) * u1
            p2 = dust_position - size_dustsheet[0]/np.sqrt(2) * u1
            p3 = dust_position + size_dustsheet[1]/np.sqrt(2) * u2
            p4 = dust_position - size_dustsheet[1]/np.sqrt(2) * u2
            ps = np.array([p1, p2, p3, p4])
            # define the x and y limits given the points just calcualted
            self.x_min = np.min(ps, axis=0)[0]
            self.x_max = np.max(ps, axis=0)[0]
            self.y_min = np.min(ps, axis=0)[1]
            self.y_max = np.max(ps, axis=0)[1]
            self.z_min = np.min(ps, axis=0)[2]
            self.z_max = np.max(ps, axis=0)[2]
    



        