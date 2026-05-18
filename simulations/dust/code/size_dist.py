import numpy as np

from astropy import units as u
from scipy.special import erf

import sys
# sys.path.append(r"C:\\Users\\tac19\\OneDrive\\Documents\\UDEL\\Project_RA\\LE\\Simulation\\code\\dust\\code")

import simulations.dust.code.var_constants as vc
# import dust_constants as dust_constants
import simulations.dust.code.fix_constants as fc

# dc = dust_constants.DustConstants(dust_env)

# curvature function
def F(a, beta, at):
    """
    Calculates the curvature term (Eq. 6) from Weingartner & Draine (2001)


    Inputs:
    -----
    a: float
    Value of the grain size in cm

    beta: float
    Constant, taken from Table 1 of Weingartner & Draine (2001)

    at: float
    Value of the upper limit for size, in cm


    Outputs:
    -----
    1 + (beta * a / at): float
    Value of F if beta is greater than or equal to zero

    (1 - (beta * a / at))^-1: float
    Value of F if beta is less than 0
    """
    
    # determine if value of beta is > or = 0
    if beta >= 0:
        # return value of first form
        return 1 + (beta * a / at)
    # if beta < 0
    else:
        # return value of second form
        return (1 - (beta * a / at))**(-1)

# D(a), size distribution for the smallest grains
def Da_carb(a, B1, B2, dc):
    """
    Calculates the value of D(a) (Eq. 2) from Weingartner & Draine (2001)


    Inputs:
    -----
    a: float
    Value of the grain size in cm

    B1: float
    Result from Bi_carb(a01, bc1)

    B2: float
    Result from Bi_carb(a02, bc2)


    Outputs:
    -----
    da1+da2: float
    Addition of the two terms from the summation
    """

    # expression for i=1
    da1 = (B1/a)*np.exp(-0.5*(np.log(a/dc.a01)/dc.sig)**2)
    # expression for i=2
    da2 = (B2/a)*np.exp(-0.5*(np.log(a/dc.a02)/dc.sig)**2)
    # sum da1+da2 and return
    return da1+da2


# including very small grains for carbonaceous dust
# Bi for B1 and B2 to go into D(a)
def Bi_carb(a0i,bci, dc):
    """
    Calculates the value of Bi (Eq. 3) from Weingartner & Draine (2001)


    Inputs:
    -----
    a0i: float
    Value of the minimum grain size in cm, a01 or a02

    bci: float
    Value of the total C abundance per H nucleus, bc1 or bc2

    NOTE: a normalization factor (sig), the density of graphite (rho), and the mass of a carbon atom (mc) are set in a separate file, dustconst.py.


    Outputs:
    -----
    Bi1*Bi2: float
    Product of the two terms for Bi, unitless
    """

    # calculate the first term of Bi
    Bi1 = ((3 / (2*np.pi)**(3/2)) * np.exp(-4.5*dc.sig**2) / (dc.rho*a0i**3*dc.sig))
    # calculate the second term of Bi
    Bi2 = bci*dc.mc / (1 + erf((3*dc.sig/np.sqrt(2)) + (np.log(a0i/dc.a01)/(dc.sig*np.sqrt(2)))))
    # multiply the two terms and return
    return Bi1*Bi2


# distribution for carbonaceous dust given grain size
def Dist_carb(a, B1, B2, dc):
    """
    Calculates the grain size distribution for carbonaceous "graphite" grains (Eq. 4) from Weingartner & Draine (2001)


    Inputs:
    -----
    a: float
    Value of the grain size in cm

    B1: float
    Result from Bi_carb(a01, bc1)

    B2: float
    Result from Bi_carb(a02, bc2)

    NOTE: constants Cg, alphag, betag, cutoff grain size (atg), and control size (acg) are set in separate file dustconst.py.


    Outputs:
    -----
    Number of grains for a particular size (units 1/cm)
    """

    # (4) = D(a) + dist1 * [1 or dist2]
    # calculate the first term, dist1
    dist1 = ((dc.Cg/a)*(a/dc.atg)**dc.alphag * F(a, dc.betag, dc.atg))

    # determine which additional term to use based on grain size
    # if a < 3.5e-8, the function is undefined so return 0
    if a.value < 3.5e-8:
        return 0
    # if 3.5e-8 <= a < atg, return the result of Da_carb(a, B1, B2) + dist1 * 1
    elif a.value >= 3.5e-8 and a < dc.atg:
        return Da_carb(a, B1, B2, dc) + dist1
    # if a > atg, return the result of Da_carb(a, B1, B2) + (dist1*dist2)
    else:
        dist2 = np.exp(-((a - dc.atg)/dc.acg)**3)
        return Da_carb(a, B1, B2, dc) + (dist1 * dist2)


# distribution for carbonaceous dust given grain size
def Dist_sili(a, dc):
    """
    Calculates the grain size distribution for carbonaceous "graphite" grains (Eq. 4) from Weingartner & Draine (2001)


    Inputs:
    -----
    a: float
    Value of the grain size in cm

    B1: float
    Result from Bi_carb(a01, bc1)

    B2: float
    Result from Bi_carb(a02, bc2)

    NOTE: constants Cg, alphag, betag, cutoff grain size (atg), and control size (acg) are set in separate file dustconst.py.


    Outputs:
    -----
    Number of grains for a particular size (units 1/cm)
    """

    # (4) = D(a) + dist1 * [1 or dist2]
    # calculate the first term, dist1
    dist1 = ((dc.Cs/a)*(a/dc.ats)**dc.alphas * F(a, dc.betas, dc.ats))

    # determine which additional term to use based on grain size
    # if a < 3.5e-8, the function is undefined so return 0
    if a.value < 3.5e-8:
        return 0
    # if 3.5e-8 <= a < atg, return the result of Da_carb(a, B1, B2) + dist1 * 1
    elif a.value >= 3.5e-8 and a < dc.ats:
        return dist1
    # if a > atg, return the result of Da_carb(a, B1, B2) + (dist1*dist2)
    else:
        dist2 = np.exp(-((a - dc.ats)/dc.acs)**3)
        return dist1 * dist2
