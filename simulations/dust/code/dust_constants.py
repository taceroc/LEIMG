# list of constants for dust analysis
# Refs:
#   Weingartner & Draine (2001); ApJ 548,296
#   Draine et al. (2021); ApJ 917,3
#   Chakradhari et al. (2019)
#   Riess et al. (2019)
#   Draine & Hensley (2021)
#   Sugerman (2003)

import astropy.units as u
# from astropy.units import cds
# cds.enable()

class DustConstants:
    def __init__(self, dust_env='mw'):
        self.dust_env = dust_env
        self.sig = 0.4
        self.rho = 2.24*u.g/u.cm**3 #g/cm^3, density of graphite
        self.mc = 1.994e-23*u.g #g, mass of 1 carbon atom
        if self.dust_env == 'mw':
            # Weingartner & Draine (2001)
            # carbonaceous dust FOR MILKY WAY
            # constants
            self.a01 = 3.5e-8*u.cm  # cm, size of grain
            self.a02 = 30.0e-8*u.cm  # cm, size of grain
            
            self.bc = 6.0e-5  # total C abund per H in log-normal pops
            self.bc1 = 0.75*self.bc
            self.bc2 = 0.25*self.bc
            
            self.alphag = -1.54
            self.betag = -0.165
            self.atg = 0.0107e-4*u.cm #cm
            self.acg = 0.428e-4*u.cm #cm
            self.Cg = 9.99e-12
            
            # silicate dust
            # constants
            self.alphas = -2.21
            self.betas = 0.300
            self.ats = 0.164e-4*u.cm  # cm
            # ats = ats.to(u.lyr)
            self.acs = 0.1e-4*u.cm  # cm
            # acs = acs.to(u.lyr)
            self.Cs = 1.00e-13
        elif self.dust_env == 'lmc':
            #FOR LMC2
            self.a01 = 3.5e-8*u.cm  # cm, size of grain
            self.a02 = 30.0e-8*u.cm  # cm, size of grain
            
            self.bc = 0.5e-5  # total C abund per H in log-normal pops
            self.bc1 = 0.75*self.bc
            self.bc2 = 0.25*self.bc
            
            self.alphag = -2.82
            self.betag = 9.01
            self.atg = 0.392e-4*u.cm #cm
            self.acg = 0.269e-4*u.cm #cm
            self.Cg = 6.20e-17
            
            # silicate dust
            # constants
            self.alphas = -2.36
            self.betas = -0.113
            self.ats = 0.182e-4*u.cm  # cm
            self.Cs = 3.03e-14
            self.acs = 0.1e-4*u.cm  # cm

# global constants
# sig = 0.4
# rho = 2.24*u.g/u.cm**3 #g/cm^3, density of graphite
# mc = 1.994e-23*u.g #g, mass of 1 carbon atom
# Rv = 3.1 #extinction, unitless
# nH = 10/u.cm**3 #atoms/cm^3, density of hydrogen
# c = 3e10*u.cm/u.s #cm/s, speed of light



#FOR LMC2
# bc = 0.5e-5  # total C abund per H in log-normal pops
# bc1 = 0.75*bc
# bc2 = 0.25*bc

# alphag = -2.82
# betag = 9.01
# atg = 0.392e-4*u.cm #cm
# acg = 0.269e-4*u.cm #cm
# Cg = 6.20e-17


# # silicate dust
# # constants
# alphas = -2.36
# betas = -0.113
# ats = 0.182e-4*u.cm  # cm
# Cs = 3.03e-14
# acs = 0.1e-4*u.cm  # cm
