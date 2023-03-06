import numpy as np
from itertools import combinations

# convert deg and hour to rad
def raddeg(d,m,s):
    return (d + m/60 + s/3600) * np.pi/180

def radhr(hr,m,s):
    return (hr + m/60 + s/3600) * np.pi/12

# the position of source Gamma Velorum 2
ra = radhr(8,9,31.95013)
dec = raddeg(-47,20,11.7108)

# the position of baseline (HESS Observatory)
lat = raddeg(-23,16,17)
lon = raddeg(16,30,00)

# the radius of aperture of telescope
rad = 6

# the effective area of each telescope
A = 100   # m^2

# Telescope position (only corner one have been taken)
#comb = combinations([-0.16-85.04j, 85.07-0.37j, 0.24+85.04j, -85.04+0.28j],2)

# best combination
base = (85.07-0.37j) - (-0.16-85.04j)                     # baseline (E,N)


# three dimensional baseline 
x = round(np.real(base),2)
y = round(np.imag(base),2)
z = 1e-6
