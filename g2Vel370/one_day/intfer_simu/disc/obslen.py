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

# wavelength of observation 
lam = 370e-9                        # disk observation (continuous spectrum)
lam1 = 465e-9                       # emission region observation (line spectrum)

# resolution time
delt = 5e-9

# observation day
day = 2460060.3                    # 25 April 2023, 7 pm


# Number of observation night (observation length)
Night = 1 

# observation duration
Nobs = 5/12                       # each for 10 hour of night

# average time of observation
M = int(Nobs * 48 + 1)                   #  dots 30 minutes apart

julday = []
for i in range(0, Night, 1):
    jday = np.linspace(day + i, day + i + Nobs, M)
    julday.append(jday)

jd = np.array(julday).flatten()     # total observation time according to Julian time

# Telescope position (only corner one have been taken)
#comb = combinations([-0.16-85.04j, 85.07-0.37j, 0.24+85.04j, -85.04+0.28j],2)

# best combination
#comb = combinations([-0.16-85.04j, 85.07-0.37j],2)      
#comb = combinations([-0.16-85.04j, -85.04+0.28j],2)

base = (-0.16-85.04j) - (85.07-0.37j)                     # baseline

