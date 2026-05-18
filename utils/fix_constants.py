# list of constants for dust analysis
# Refs:
#   Weingartner & Draine (2001); ApJ 548,296
#   Draine et al. (2021); ApJ 917,3
#   Sugerman (2003)

# NOTE: All constants with units are in year and light years
dtoy = 0.00273973 # 1 day = 0.00273973 y
pctoly = 3.26156 # 1pc = 3.26156 light-year
pctom = 3.086e+16 # 1pc = 3.086e16 meter
lsuntol = 3.9e36 # L sun  # watts = kg m ^2 / s^3
ytos = 3.154e+7 # 1y = 3.154e+7s
c = 1 # light speed in ly / y
n_H = 300000 #* (100 ** 3) #cm-3
#n_H = n_H * ( pctom ** 3 ) / ( pctoly ** 3 ) # ly-3
sigma = 5e-22 / (100 ** 2) #m2 #RR paper
sigma = ( sigma / ( (pctom ** 2) ) ) * ( pctoly ** 2 )
albedo = 0.6

Msun = 4.75 # absolute magnitde sun