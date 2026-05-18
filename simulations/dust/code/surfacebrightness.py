"""
Functions needed to calculate the surface brightness of scattered light (Eq. 7) from Sugerman (2003); ApJ 126,1939.
"""

# __all__ = [Flamb,
#            deltap,
#            Bscat]
__author__ = 'Charlotte M. Wood'
__version__ = '0.1'


import numpy as np
import astropy.units as u

import dustconst as c
import grainsizedist as gsd
import interpdust as id
import scatintegral as scaint


# calculate F(lambda)
def Flamb(Ftmax):
    """
    Calculates the time-integrated monochromatic flux, F(lambda). See paragraph after Eq. 15 in Sugerman (2003); ApJ 126,1939.


    Inputs:
    -----
    Ftmax: float
    Value of the flux at a given wavelength for the SN at peak, in erg/s/cm^2

    NOTE: 0.5 * tau is substituted for t2, since tau = 2*t2. Tau is the time from +2mag pre-peak to +2mag post-peak and is set in a separate file, dustconst.py


    Outputs:
    -----
    Flamb: float
    Value of the time-integrated monochromatic flux, in erg/cm^2
    """

    # F(lambda) = 1.25 F(lambda, tmax) * t2
    # F(lambda) = 1.25 Ftmax * 0.5 tau
    Flamb = 1.25*Ftmax*0.5*c.tau.to(u.s)
    return Flamb


# calculate delta p
def deltap(deltaz):
    """
    Calculate the width of the observed light echo, Delta p. See Eq. 11 and Fig. 2 of Sugerman (2003); ; ApJ 126,1939.


    Inputs:
    -----
    deltaz: float
    Value of the width of the dust sheet causing the light echo in cm. See Fig. 2 of Sugerman (2003); ApJ 126,1939

    NOTE: The following values are set in a separate file, dustconst.py
        - pecho: the physical extent of the light echo in cm (observed angular size * distance)
        - t: time of observation, in days since peak
        - c: the speed of light in cm/s
        - tau: time between +2mag pre-peak and +2mag post-peak
    
    
    Outputs:
    -----
    deltap: float
    Value of the width of the observed light echo in cm
    """

    deltap = np.sqrt((((c.pecho.to(u.cm)/(2*c.t.to(u.s))) + ((c.c**2*c.t.to(u.s))/(2*c.pecho.to(u.cm))))**2 * c.tau.to(u.s)**2) + ((c.c*c.t.to(u.s)/c.pecho.to(u.cm))**2 * deltaz**2))
    return deltap


# calculate the surface brightness
def Bscat(Ftmax, deltaz, S, radcm):
    """
    Calculates the surface brightness of scattered light at a given wavelength and time of observation, Bsc(lambda, t, phi). See Eq. 7 of Sugerman (2003); ApJ 126,1939.


    Inputs:
    -----
    Ftmax: float
    Value of the flux at a given wavelength for the SN at peak, in erg/s/cm^2

    deltaz: float
    Value of the width of the dust sheet causing the light echo in cm. See Fig. 2 of Sugerman (2003); ApJ 126,1939

    S: float
    Value of the integrated scattering function across all sizes for a given wavelength. Result of integrating dS from scatintegral.py. 

    radcm: float
    r, value of the strait-line distance between the SN and the scattering dust in units of cm. See Fig. 1 of Sugerman (2003); ApJ 126,1939.

    NOTE: the number density of hydrogen atoms (nH), the speed of light in cm/s (c), and the physical extent of the light echo in cm (pecho; observed angular size * distance) are set in a separate file, dustconst.py.


    Outputs:
    -----
    Bscat: float
    Value of the surface brightness of scattered light at a given wavelength and time of observation, in units of flux/area (erg/s/cm^2/cm^2)
    """

    # Bsc(lambda, t, phi) = F(lambda) * nH(r) * [c deltaz / 4 pi r p deltap] * S(lambda,  mu)
    Bscat = Flamb(Ftmax) * c.nH * c.c * deltaz * S / (4*np.pi*radcm*c.pecho.to(u.cm)*deltap(deltaz))
    return Bscat