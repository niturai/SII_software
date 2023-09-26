import numpy as np
import math
from scipy.optimize import leastsq
from julian import dst, den

# path to read file
drk = '../data/orbit/'

# the general constant value
c = 299792458.0                                  # in m/s
h = 6.62607015e-34                               # J.s
k = 1.380649e-23                                 # J/K

# the constant parameter of a star
T_a = 35000                                      # Temperature of source A
T_b = 57000                                      # Temperature of disk
T1 = 13880                                       # Temperature of atmosphere of WR
R_a = 39.4503                                    # Radius of source A in light seconds
R_b = 13.9236                                    # Radius of disk
D = 3.458e10                                     # Distance from earth in light seconds

# position of orbit in sky with time
dAB = np.load(drk + "dAB.npy")
Xth = []
Yth = []
for i in range(0, len(dst)):
    Xth.append(dAB[0, dst[i] : den[i]])
    Yth.append(dAB[1, dst[i] : den[i]])

X = np.concatenate(Xth)                          # X-position of orbit according to observation time each day
Y = np.concatenate(Yth)                          # Y-position of orbit according to observation time each day

# define the radius of emission region, considering the first position of binary system
def emis(M1,M2):
    A = np.sqrt(X[0]**2 + Y[0]**2)               # orbital separation of binary
    q = M1/M2                                    # mass ratio (M1 having the roche lobe)
    a1 = 0.49*q**(2/3)
    a2 = 0.6*q**(2/3)
    a2 += math.log(1+q**(1/3))
    r1 = A*a1/a2
    return r1

R1 = emis(9,28.5) - R_b                          # Radius of WR (disk + atmosphere)
print('Radius of roche lobe from surface',R1)

