# list of constants for dust analysis
# Refs:
#   Weingartner & Draine (2001); ApJ 548,296
#   Draine et al. (2021); ApJ 917,3
#   Chakradhari et al. (2019)
#   Riess et al. (2019)
#   Draine & Hensley (2021)
#   Sugerman (2003)

# NOTE: All constants with units are in cgs units
import astropy.units as u
from astropy.units import cds
cds.enable()

# global constants
sig = 0.4
rho = 2.24*u.g/u.cm**3 #g/cm^3, density of graphite
mc = 1.994e-23*u.g #g, mass of 1 carbon atom
Rv = 3.1 #extinction, unitless
nH = 10/u.cm**3 #atoms/cm^3, density of hydrogen
c = 3e10*u.cm/u.s #cm/s, speed of light

# 09ig constants
tau = 37.25*u.d #days between +2mag before peak to +2mag after peak
tpeak = 0.0*u.d #days since 55080.14 MJD, t~ in Sugerman (2003)
t = 1446.83*u.d #days since 55080.14 MJD (56526.97-55080.14)
dist = 103.0e6*u.lyr #ly, distance to SN from earth
disterr = 4.0e6*u.lyr #ly, error on distance to SN
z = 130.0*u.lyr #ly, distance between SN and echo dust
zerr = 6.0*u.lyr #ly, error on distance to dust
pecho = 32.3*u.lyr #ly, physical extent of echo (angular size * distance)
pechoerr = 0.7*u.lyr #ly, error on extent of echo

# Weingartner & Draine (2001)
# carbonaceous dust
# constants
a01 = 3.5e-8*u.cm  # cm, size of grain
a02 = 30.0e-8*u.cm  # cm, size of grain
bc = 6.0e-5  # total C abund per H in log-normal pops
bc1 = 0.75*bc
bc2 = 0.25*bc
alphag = -1.54
betag = -0.165
atg = 0.0107e-4*u.cm #cm
acg = 0.428e-4*u.cm #cm
Cg = 9.99e-12

# silicate dust
# constants
alphas = -2.21
betas = 0.300
ats = 0.164e-4*u.cm  # cm
acs = 0.1e-4*u.cm  # cm
Cs = 1.00e-13

# Draine et al. (2021) "astrodust" & PAH
# constants
p = -3.3
aminAd = 4.5e-8*u.cm #cm, min size of astrodust particles
amaxAd = 0.4e-4*u.cm #cm, max size of astrodust particles
Vad = 5.34e-27*u.cm**3 #cm^3/H, volume of dust
poro = 0.2 #porosity of grains
aminPAH = 4.0e-8*u.cm #cm, min size of PAH particles
amaxPAH = 1e-6*u.cm #cm, max size of PAH particles
a01pah = 4.0e-8*u.cm #cm
a02pah = 30e-8*u.cm #cm
B1pah = 6.134e-7
B2pah = 3.113e-10
Vpah1 = 3.0e-28*u.cm**3 #cm^3/H, volume of dust
Vpah2 = 0.7e-28*u.cm**3 #cm^3/H, volume of dust

# Hensley & Draine (2023) astrodust + PAH
# constants
B1 = 7.52e-7 #1/H
B2 = 8.09e-10 #1/H
BAd = 3.31e-10 #1/H
a0Ad = 63.8e-8 * u.cm #cm
sigAd = 0.353
aminAd23 = 4.5e-8*u.cm #cm, min size for astrodust particles
amaxAd23 = 5.0e-4*u.cm #cm, max size for astrodust particles
A0 = 2.97e-5 #1/H
Ai = [-3.40, -0.807, 0.157, 7.96e-3, -1.68e-3]
A1 = -3.40
A2 = -0.807
A3 = 0.157
A4 = 7.96e-3
A5 = -1.68e-3