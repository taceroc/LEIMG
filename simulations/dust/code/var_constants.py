# list of constants for dust analysis
# Refs:
#   Weingartner & Draine (2001); ApJ 548,296
#   Draine et al. (2021); ApJ 917,3
#   Sugerman (2003)

# Variable constants

import simulations.dust.code.fix_constants as cf
import numpy as np

Deltat = 180 # days #time of LE detection after source detection/peak detection
Deltat_y = np.array(Deltat) * cf.dtoy
# https://academic.oup.com/mnras/article/480/2/1466/5065048
dt0 = 0.05 * cf.dtoy #eta carinae # duration light from source in years
dz0 = 0.02 * cf.pctoly # in ly # thickness of dust sheet
ct = cf.c * Deltat_y


z0 = 0.1 # pc # distance source-dust
z0ly = np.array(z0) * cf.pctoly
alpha = 18 # inclination plane
r0 = 2 #pc # radii of the dust sphere
r0ly = np.array(r0) * cf.pctoly


dkpc = 5 #kpc #distance source-earth
d = 19000 #dkpc * 1000 * cf.pctoly
# L = 15,000 x 3.9e26 # watts = kg m ^2 / s^3
L = 15e3
L = (L * cf.lsuntol) * (cf.ytos ** 3) * (cf.pctoly ** 2) / (cf.pctom ** 2)
m_peak = 6.75