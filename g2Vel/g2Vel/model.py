import numpy as np
import math
import scipy.special as sp
from starpara import c, h, k, T_a, T_b, T1, R_a, R_b, R1, D, t, zs

# orbit fitting coefficients for X and Y
x0, x1, x2, x3 = zs[0], zs[2], zs[4], zs[6]
y0, y1, y2, y3 = zs[1], zs[3], zs[5], zs[7]


# photon flux from single source
def pho_spec(R, T, lam):
    a1 = np.pi * R**2
    a1 /= (D*lam)**2
    a2 = h*c/(lam*k*(T + 2.725))
    a2 = np.clip(a2,0,10)
    I = a1 / (np.exp(a2) - 1)
    return I

# define the model of signal
def hbt(x, y, t, lam, lam1, R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3):           # lam is for continuum and lam1 is for line spectrum

    I_a = pho_spec(R_a, T_a, lam)
    I_b = pho_spec(R_b, T_b, lam)
 
    if lam1 == lam:
       I_1 = pho_spec(R1, T1, lam1)
       I_b1 = pho_spec(R_b, T1, lam1)
    else:
       I_1 = I_b1 = 0

    r = np.sqrt(x**2 + y**2)
    v = 2*np.pi/(lam*D)
    
    b1 = v*r*R_a
    b2 = v*r*R_b
    b3 = v*r*R1
    b4 = v*r*R_b

    X = x0 + x1*t + x2*t**2 + x3*t**3
    Y = y0 + y1*t + y2*t**2 + y3*t**3

    V_a = 2 * I_a * sp.jv(1, b1)/b1

    V_2 = I_b * sp.jv(1, b2)/b2
    V_3 = I_1 * sp.jv(1, b3)/b3
    V_4 = I_b1 * sp.jv(1, b4)/b4

    V_b = 2 * (V_2 + V_3 - V_4)
    
    V_ab = V_a * V_b * np.cos(v * (x*X + y*Y))

    V = (V_a**2 + V_b**2 + 2 * V_ab)
    V /= (I_a + I_b + I_1 - I_b1)**2
    
    return V


