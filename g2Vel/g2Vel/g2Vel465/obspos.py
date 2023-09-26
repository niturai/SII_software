import numpy as np
from itertools import combinations

# convert deg and hour to rad
def raddeg(d,m,s):
    return (d + m/60 + s/3600) * np.pi/180

def radhr(hr,m,s):
    return (hr + m/60 + s/3600) * np.pi/12

# the position of Regor
ra = radhr(8,9,31.95013)
dec = raddeg(-47,20,11.7108)

# the position of HESS Observatory
lat = raddeg(-23,16,17)
lon = raddeg(16,30,00)

# the radius of aperture of telescope
rad = 6

# the effective area (in m^2) of each telescope
A = 100  

'''
# Telescope position (C1, C2, C3, C4)
comb = combinations([85.04 - 0.16j, 0.37 + 85.07j, -85.04 + 0.24j, -0.28 - 85.04j], 2)          # (X, Y) = (E, N)

# Number of baselines for observations
base = []
for i in list(comb):
    base.append(i[1] - i[0])
'''

# selected telescope
base = (-0.28 - 85.04j) - (-85.04 + 0.24j)   

# three dimensional baseline 
x = np.real(base)
y = np.imag(base)
z = 1e-6

