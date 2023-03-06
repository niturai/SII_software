import numpy as np
from starpara import R1, R_b, X, Y

# flux on carbon line
Phi = 2e-7              # photons m^-2 s^-1 Hz^-1
d = 1.03668231976e+19   # in meter

# luminosity
L = 4 * np.pi * d**2 * Phi

print(L)

# radius of the WR
R = (R1 + R_b) * 2.998e+8   # in meter

# radius of the WR + o system
R = np.sqrt(X[0]**2 + Y[0]**2) * 2.998e+8       # in meter

# Stefan-Boltzmann constatnt
sigma = 5.67037442e-8 

# temperature
T_eff = (L/(4 * np.pi * R**2 * sigma))**(1/4)

print(T_eff)
