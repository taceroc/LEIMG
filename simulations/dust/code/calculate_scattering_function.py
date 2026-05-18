import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy import units as u
from scipy.special import erf
from scipy import integrate

import sys
import os

### credits to https://github.com/cwood12/lightecho-dust
import simulations.dust.code.var_constants as vc
import simulations.dust.code.dust_constants as dust_constants
import simulations.dust.code.fix_constants as fc
import simulations.dust.code.scattering_function as sf
import simulations.dust.code.size_dist as sd


def calculate_scattering_function_values(wave, sizeg, waveg, Qcarb, gcarb,
                                        Qsili, gsili, dust_env):
    dc = dust_constants.DustConstants(dust_env)
    B1 = sd.Bi_carb(dc.a01, dc.bc1, dc)
    # calculate B2
    B2 = sd.Bi_carb(dc.a02, dc.bc2, dc)
    carbon_distribution = [sd.Dist_carb(idx, B1, B2, dc).value for idx in 1e-4*sizeg*u.cm] #in cm
    silicone_distribution = [sd.Dist_sili(idx, dc).value for idx in 1e-4*sizeg*u.cm] #in cm
    print("carb values")
    Qc_scs, gc_s, sizes, carbon_distribution = sf.dS_pre_interpolation(sizeg, waveg, wave, Qcarb, gcarb, carbon_distribution)
    print("silicate values")
    Qs_scs, gs_s, sizes, silicone_distribution = sf.dS_pre_interpolation(sizeg, waveg, wave, Qsili, gsili, silicone_distribution)
    # print("carb values")
    # Qc_scs, gc_s, sizes, carbon_distribution = sf.dS_pre(sizeg, waveg, wave, Qcarb, gcarb, carbon_distribution)
    # print("silicate values")

    # Qs_scs, gs_s, sizes, silicone_distribution = sf.dS_pre(sizeg, waveg, wave, Qsili, gsili, silicone_distribution)


    return  sizes, Qc_scs, gc_s, carbon_distribution, Qs_scs, gs_s, silicone_distribution

def calculate_scattering_function(mu, sizes, Qc_scs, gc_s, carbon_distribution,
                                  Qs_scs, gs_s, silicone_distribution, composition='both'):

    # Qscs, gs, sizes, carbon_distribution = calculate_scattering_function_values(mu, sizeg, waveg, wave, Qcarb, gcarb)
    
    ds_c = sf.dS(mu, Qc_scs, gc_s, sizes, carbon_distribution)
    S_c = sf.S(ds_c, sizes)

    ds_s = sf.dS(mu, Qs_scs, gs_s, sizes, silicone_distribution)
    S_s = sf.S(ds_s, sizes)

    assert composition in ['both', 'S', 'C'], "Composition is not valid, must be 'both', 'S' or 'C'"

    if composition == 'both':
        return ds_c+ds_s, S_c+S_s
    elif composition == 'S':
        return ds_s, S_s
    elif composition == 'C':
        return ds_c, S_c

def load_data():
    # path_dustdata = '/content/drive/MyDrive/LE2023/dust/data/'
    # path_dustdata = r"C:\\Users\\tac19\\OneDrive\\Documents\\UDEL\\Project_RA\\LE\\Simulation\\code\\dust\\data"
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) 
    # print(ROOT_DIR)
    BASE_DIR = os.path.join( os.path.dirname( __file__ ), '..' )
    # print(BASE_DIR)
    # PATH_TO_DUST_IMG = os.path.join(ROOT_DIR, 'dust/cube/Data') 
    path_dustdata = os.path.join(BASE_DIR, 'data')
    # print(path_dustdata)
    # path_dustdata = r"C:\\Users\\tac19\\OneDrive\\Documents\\UDEL\\Project_RA\\LE\\lightecho_modeling_oop\\OOP\\dust\\data" 
    # pull out available wavelengths for g values, convert to cm from um, and take the 
    waveg = np.loadtxt(path_dustdata+r'/dustmodels_WD01/LD93_wave.dat', unpack=True) #micronm
    # pull out available sizes for the g values, convert to cm from um, and take the log
    sizeg = np.loadtxt(path_dustdata+r'/dustmodels_WD01/LD93_aeff.dat', unpack=True) #micron


    # carbonaceous dust
    carbonQ = path_dustdata+r'/dustmodels_WD01/Gra_81.dat'
    Qcarb_sca = np.loadtxt(carbonQ, usecols=(2), unpack=True)
    Qcarb_abs = np.loadtxt(carbonQ, usecols=(1), unpack=True)
    Qcarb = Qcarb_sca / (Qcarb_sca + Qcarb_abs)
    # g (degree of forward scattering) values
    gcarb = np.loadtxt(carbonQ, usecols=(3), unpack=True)

    # silicate dust
    siliconQ = path_dustdata+r'/dustmodels_WD01/suvSil_81.dat'
    # Qsil = np.loadtxt(siliconQ, usecols=(2), unpack=True)
    Qsili_sca = np.loadtxt(siliconQ, usecols=(2), unpack=True)
    Qsili_abs = np.loadtxt(siliconQ, usecols=(1), unpack=True)
    Qsili = Qsili_sca / (Qsili_sca + Qsili_abs)
    # g (degree of forward scattering) values
    gsili = np.loadtxt(siliconQ, usecols=(3), unpack=True)


    return sizeg, waveg, Qcarb, gcarb, Qsili, gsili


def main(mu, sizes, Qc_scs, gc_s, carbon_distribution, Qs_scs, gs_s, silicone_distribution, composition):
    # sizeg, waveg, Qcarb, gcarb = load_data

    ds, S = calculate_scattering_function(mu, sizes, Qc_scs, gc_s, carbon_distribution,
                                  Qs_scs, gs_s, silicone_distribution, composition)

    return ds, S