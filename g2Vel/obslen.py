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
lam = 6e-7

# resolution time
delt = 5e-9

# observation day
day = 2460000                  # 24 feb 2023

# sqrt(channels x nights)
chnl = 1

# observation duration
Nobs = 1/2                       # each for 12 hour of night

# Number of observation night (observation length)
Night = 1 

# average time of observation
M = int(Nobs * 86400 + 1)                   #  dots 1 sec apart

julday = []
for i in range(0, Night, 1):
    jday = np.linspace(day + i, day + i + Nobs, M)
    julday.append(jday)

jd = np.array(julday).flatten()     # total observation time according to Julian time

# Telescope position
comb = combinations([-60+60j, 60+60j, -60-60j, 60-60j],2)

# Number of baselines for observations
base = []
for i in list(comb):
    base.append(i[1] - i[0])

base = base


# the circular aperture of each telescope
def circ(x,y,rad):                            # rad is the radius of telescope
    f = 0*x
    f[x*x + y*y < rad**2] = 1
    return f


