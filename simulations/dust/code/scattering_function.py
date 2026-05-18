import numpy as np

from astropy import units as u
from scipy.special import erf
from scipy import integrate, interpolate

import sys
sys.path.append(r"C:\\Users\\tac19\\OneDrive\\Documents\\UDEL\\Project_RA\\LE\\Simulation\\code\\dust\\code")

import simulations.dust.code.var_constants as vc
import simulations.dust.code.dust_constants as dc
import simulations.dust.code.fix_constants as fc


# extract Qsc(lambda, a) value from table
def extract_Qsc_g(wave, a, sizeg, waveg, valuesq):
    """
    Extracts and returns the value of Qsc from interpolated data table from Draine et al. (2021); ApJ 917,3. Qsc describes the grain scattering efficiency.

    
    Inputs:
    -----
    wave: float
    A specified value for the wavelength in cm, on which Qsc depends

    a: float
    A specified value for the size of the grain in cm, on which Qsc depends

    pointsq: array-like
    An array of (wavelength, size) points

    valuesq: array-like
    An array of Qsc values corresponding to each (wavelength, size) point

    
    Outputs:
    -----
    valuesq[idx]: float
    The corresponding value of Qsc for the specified (wave, a) point
    """

    # a = 1.259E-03
    # w = 6.310E+02
    # iterate over the array of points to search for matching (wave, a) point
    for idx in range(0,len(sizeg)):
        for idy in range(0,len(waveg)):
            # if the first value of the current point matches wave and if the second value matches a
            # return the corresponding Qsc value
            if sizeg[idx] == a and waveg[idy] == wave:
                ids = len(waveg) * idx + idy
                return valuesq[ids]
            # if one or both do not match, move to the next iteration
            else:
                continue


# calculate the phase function Phi(mu, lambda, a)
def Phi(mu, g):
    """
    Calculates the scattering phase function, Phi(mu, lambda, a). See Eq. 3 of Sugerman (2003); ApJ 126,1939.

    
    Inputs:
    -----
    mu: float
    Value of cos(scattering angle) for a given epoch of observation

    g: float
    Value of the degree of forward scattering for a given size and wavelength (from extract_g)


    Outputs:
    -----
    phi: float
    Value of the scattering phase function for the given scattering angle, wavelength, and size
    """

    phi = (1 - g**2) / ((1 + g**2 - 2*g*mu)**(3/2))
    return phi


# def dS(mu, sizeg, waveg, wave, Qcarb, gcarb, distribution):
#     ds = []
#     # wave = 6.310E+02
#     for idy in range(0,len(waveg)):
#         if waveg[idy] == wave:
#             indy = idy
#     # iterate over the array of points to search for matching (wave, a) point
#     # for m in mu:
#     for idx in range(0,len(sizeg)):
#         # if the first value of the current point matches wave and if the second value matches a
#         # return the corresponding Qsc value
#         ids = len(waveg) * idx + indy
#         Qsc = Qcarb[ids]
#         g = gcarb[ids]
#         phi = Phi(mu, g)
#         #carbon_distribtion
#         f = distribution[idx]
#         ds.append((Qsc * np.pi*(1e-4*sizeg[idx])**2 * phi * f)) #in cm

#     return ds

def dS_pre(sizeg, waveg, wave, Qcarb, gcarb, distribution):
    Qscs, sizes, dist, gs = [], [], [], []
    # wave = 6.310E+02
    indy = 0
    for idy in range(0,len(waveg)):
        if waveg[idy] == wave:
            indy = idy
    # iterate over the array of points to search for matching (wave, a) point
    # for m in mu:
    for idx in range(0,len(sizeg)):
        # if the first value of the current point matches wave and if the second value matches a
        # return the corresponding Qsc value
        ids = len(waveg) * idx + indy
        Qsc = Qcarb[ids]
        g = gcarb[ids]
        # phi = Phi(mu, g)
        #carbon_distribtion
        # f = distribution[idx]
        Qscs.append(Qsc)
        sizes.append(sizeg[idx])
        dist.append(distribution[idx])
        gs.append(g)
        # ds.append((Qsc * np.pi*(1e-4*sizeg[idx])**2 * phi * f)) #in cm
        # print(idx)
        # print("VALUESS DUSTT")
        # print(waveg[indy] ,Qsc, ids, sizeg[idx])
    return Qscs, gs, sizes, dist


def dS_pre_interpolation(sizeg, waveg, wave, Qcarb, gcarb, distribution):
    Qscs, sizes, dist, gs = [], [], [], []
    # wave = 6.310E+02
    # for idy in range(0,len(waveg)):
    #     if waveg[idy] == wave:
    #         indy = idy
    # iterate over the array of points to search for matching (wave, a) point
    # for m in mu:
    for idx in range(0,len(sizeg)):
        # if the first value of the current point matches wave and if the second value matches a
        # return the corresponding Qsc value
        ids = len(waveg) * idx #+ len(waveg)
        Qsc_all = []
        g_all = []
        for i in range(ids, ids+len(waveg))[::-1]:
            Qsc_all.append(Qcarb[i])
            g_all.append(gcarb[i])
        spl = interpolate.make_interp_spline(np.array(waveg)[::-1], np.c_[np.array(Qsc_all), np.array(g_all)])
        Qsc, g = spl(wave).T
        Qscs.append(Qsc)
        sizes.append(sizeg[idx])
        dist.append(distribution[idx])
        gs.append(g)
        # ds.append((Qsc * np.pi*(1e-4*sizeg[idx])**2 * phi * f)) #in cm
        # print("VALUESS DUSTT")
        # print(Qsc, sizeg[idx])
    return Qscs, gs, sizes, dist

def dS(mu, Qsc, g, sizes, f):
    ds = []
    phi = Phi(mu, np.array(g))
    # wave = 6.310E+02

    ds.append((np.array(Qsc) * np.pi*(1e-4*np.array(sizes))**2 * np.array(phi) * np.array(f))) #in cm

    return ds

def S(ds, sizeg):
    S = []
    # for i in range(len(mu)):
    # S.append(integrate.simpson(ds, 1e-4*np.array(sizeg))) # in cm
    S.append(integrate.trapezoid(ds, 1e-4*np.array(sizeg))) # in cm


    return S