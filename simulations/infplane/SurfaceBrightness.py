import numpy as np
import sys
from utils import fix_constants as fc
from simulations.dust.code import calculate_scattering_function as csf
import logging
import pandas as pd
logger = logging.getLogger(__name__)

class SurfaceBrightness:
    def __init__(self, wavel, source, LE, composition_s_c='both', dust_env='mw', lc=None):
        """
        Calculate the surface brightness at a position r = (x_inter, y_inter, z_inter):
        Sugermann 2003 equation 7:
            SB(lambda, t) = F(lambda)nH(r) * (c dz0 / (4 pi r rhodrho) )* S(lambda, mu)
            S(lambda, mu) = \int Q(lamdda, a) sigma Phi(mu, lambda, a) f(a) da
            lambda: given wavelength in micrometer [1/lenght]
            dz0: dust thickness [lenght]
            r: position dust [lenght]
            rhodrho: x-y of LE [lenght^2]
            mu: cos theta, theta: scattering angle
            Q: albedo
            sigma: cross section [lenght^2]
            Phi: scattering function
            f(a): dust distribution [1/lenght]
            S: scattering integral [lenght^2]

            Return units [mass / (time^2 lenght^3)] >> [flux / lenght^3]
        """
        self.wavel = wavel
        self.Fl = source.Flmax
        self.dt0 = source.dt0
        self.d = source.d
        self.ct = LE.ct
        self.dz0 = LE.dz0
        self.r_le2 = LE.calculate_rle2()
        self.rhos = 0
        self.surface = 0
        self.cossigma = 0
        self.lc = lc
        self.sb_true_matrix = 0
        self.composition = composition_s_c
        self.dust_env = dust_env

    def define_bandpass_rubin(self):
#         {'band_name': ['lsstu', 'lsstg', 'lsstr', 'lssti', 'lsstz', 'lssty'],
        lsstu = [np.array([3104.9999999999995, 4085.9999999999995]) * 1e-4, 'u', 24.5, 10**((-48.6 - 24.5) / 2.5)]
        lsstg = [np.array([3865.9999999999995, 5669.999999999999]) * 1e-4, 'g', 25, 10**((-48.6 - 25) / 2.5)]
        lsstr = [np.array([5369.999999999999, 7059.999999999999]) * 1e-4, 'r', 25, 10**((-48.6 - 25) / 2.5)]
        lssti = [np.array([6759.999999999999, 8329.999999999998]) * 1e-4, 'i', 24, 10**((-48.6 - 24) / 2.5)]
        lsstz = [np.array([8029.999999999998, 9384.999999999998]) * 1e-4, 'z', 23.5, 10**((-48.6 - 23.5) / 2.5)]
        lssty = [np.array([9083.999999999998, 10944.999999999998]) * 1e-4, 'y', 23.5, 10**((-48.6 - 23.5) / 2.5)]

        for ixi, bandpasses in enumerate([lsstu, lsstg, lsstr, lssti, lsstz, lssty]):
            if bandpasses[0][0] <= self.wavel <= bandpasses[0][1]:
                # band = bandpasses[1]
                self.band_pass_index = ixi


    def rhos_half(self):
        """
            Calculate the thickness of the visible light echo, the rho coordiante, Sugermann 2003. Eq 11
            Convolution of the thickness due to dust thickness and duration of pulse from source
            Arguments:
                Values of z: intersection paraboloid+dust
            
            Return:
                rhodrho: Sugermann 2003 Eq 7, rho = sqrt(x**2 + y**2)
                rhos = sqrt(x**2 + y**2)
                half_obs_thickness = thickness of LE

        """
        self.rhos = np.sqrt(2 * self.z_inter_values * self.ct + (self.ct) ** 2)
        # self.rhos = np.sqrt(self.x_inter_values**2+self.y_inter_values**2+self.z_inter_values**2) * np.sin(np.deg2rad(35))

        self.half_obs_thickness = (
            np.sqrt((self.ct / self.rhos) ** 2 * self.dz0**2
                + ((self.rhos * fc.c / (2 * self.ct)) + (fc.c * self.ct / (2 * self.rhos))) ** 2 * self.dt0**2)/ 2
        )
        # self.half_obs_thickness = np.sqrt( (self.ct / self.rhos) ** 2 * self.dz0 ** 2 )
        self.rhodrho = self.rhos * self.half_obs_thickness
        return self.rhodrho, self.rhos, self.half_obs_thickness
    
    def light_curve_integral(self, tilde=0):
        """
            Calculate integral below Sugermann 2003 equation 5.
            Integral of the light curve F(lambda) = \int F(lamnda, t_tilde) dt_tilde
            t_tilde: time at the dust position (state of light curve when light interact with dust at a tiem t_tilde)
            t_tilde: between eq 4 and 5 from Sugermann 2003 (also in Xu & Crotts, 1994ApJ 435)

            Arguments:
                Values of z: intersection paraboloid+dust
                light curve: {'time': values in days centered at peak, 'mag': values in mag}
            
            Return:
                Flux at wavelenght lambda
        """
        self.Fl = 0 
        lc_integral_time = self.lc['time'][self.lc['time'] <= tilde] * 3.154E+7 #in seconds
        # lc_integral_time = self.lc['time'] * 3.154E+7 #in seconds

        self.diffetime = abs(lc_integral_time.min() - lc_integral_time.max())

        lc_integral_mag = self.lc['mag'][self.lc['time'] == tilde].values[0] - 5
        # lc_integral_mag = self.lc['mag']
        flux_upto = 10**((-48.6 - lc_integral_mag) / 2.5)
        self.Fl = flux_upto * 3e10 / ((self.wavel*1e-4))**2 * (1e-5) #(3 × 10¹⁰ cm/s / (5.5 × 10⁻⁵ cm)²) × (10⁻⁵ cm)
        logger.info("FLUX INTEGRAL IN TIME %s", self.Fl)
  
    def load_dust_values(self):
        sizeg, waveg, Qcarb, gcarb, Qsili, gsili= csf.load_data()
        # Calculate the scattering integral and the surface brightness
        (sizes, Qc_scs, gc_s, carbon_distribution, 
         Qs_scs, gs_s, silicone_distribution) = csf.calculate_scattering_function_values(
            self.wavel, sizeg, waveg, Qcarb, gcarb, Qsili, gsili, self.dust_env
        )

        return (sizes, Qc_scs, gc_s, carbon_distribution, 
                Qs_scs, gs_s, silicone_distribution)
    

    def determine_flux_time_loop(self, tilde=0):
        self.Ir = 0 #np.ones(len(r))
        if self.lc.shape[0] == None:
            # self.rhos_half()
            Fl = self.Fl #* (fc.ytos**3)  # kg,ly,y
            self.Ir = self.Ir * Fl * fc.n_H * fc.c#*1.25 * 0.5 * self.dt0  #* fc.c
        else:
            self.light_curve_integral(tilde)
            Fl = np.array(self.Fl) #* (fc.ytos**2)  # kg,ly,y
            self.Ir = Fl * fc.n_H 

class SurfaceBrightnessAnalytical(SurfaceBrightness):
    def __init__(self, wavel, source, LE, xyz_intersection, composition_s_c='both', dust_env='mw'):
        super().__init__(wavel, source, LE)
        self.x_inter_values = xyz_intersection[0]
        self.y_inter_values = xyz_intersection[1]
        self.z_inter_values = xyz_intersection[2]
    
        

    def calculate_surface_brightness(self, tilde):
        """
        Arguments:
            x_inter, y_inter, z_inter: intersection paraboloid + dust in ly
            dz0: thickness dust in ly
            ct: time where the LE is observed in y

        Return
            Surface brightness in units of kg/(ly^2 y^2) ###kg/ly^3 [erg/(s cm4)]
            cos(scatter angle)
        """

        # calculate r, source-dust
        r = np.sqrt(
            self.x_inter_values**2 + self.y_inter_values**2 + self.z_inter_values**2
        )

        # Sugerman 2003 after eq 15 F(lambda) = 1.25*F(lambda, tmax)*0.5*dt0
        # F 1.08e-14 # watts / m2
        # Ir = np.ones(len(r))
        if self.lc.shape[0] == None:
            super().determine_flux_time_loop()
        else:
            super().determine_flux_time_loop(tilde)
        

        self.sb_true_matrix = np.zeros(len(r))
        rhodrho, rhos, half_obs_thickness = super().rhos_half()
        logger.info('mean of rho %s', np.mean(rhos))
        logger.info('half_obs_thickness given by dust %s', np.mean(half_obs_thickness))
  
        # dust-observer
        ll = np.sqrt(
            self.x_inter_values**2
            + self.y_inter_values**2
            + (self.z_inter_values - self.d) ** 2
        )
        # calcualte scatter angle, angle between source-dust , dust-observer
        self.cossigma = (
            self.x_inter_values**2
            + self.y_inter_values**2
            + self.z_inter_values * (self.z_inter_values - self.d)
        ) / (r * ll)
        S = np.zeros(len(r))
        (sizes, Qc_scs, gc_s, carbon_distribution, 
         Qs_scs, gs_s, silicone_distribution) = super().load_dust_values()
        # self.cossigma = np.ones_like(self.cossigma)*np.cos(np.deg2rad(35))
        for ik, rm in enumerate(self.cossigma):
            if (rm >= -1) and (rm <= 1):
                ds, Scm = csf.main(
                    rm, sizes, 
                    Qc_scs, gc_s, carbon_distribution, 
                    Qs_scs, gs_s, silicone_distribution, composition=self.composition
                )  # 1.259E+00 in um
                S[ik] = (Scm[0][0])# * fc.pctoly**2) / ((100 * fc.pctom) ** 2 )

            else:
                S[ik] = 0
  
        Inte_z = self.dz0 * 9.461e+17 # in cm

        self.sb_true_matrix = (
            self.Ir * S * (Inte_z) * (2.998e+10) / (4 * np.pi * (r * 9.461e+17)) # AB mag units erg, cm, s
        )  
      
        logger.info("FLUX AFTER SURFACE %s", np.mean(self.sb_true_matrix))#/(rhodrho* 9.461e+17**2)))
        return self.cossigma, self.sb_true_matrix#/(rhodrho* 9.461e+17**2)
