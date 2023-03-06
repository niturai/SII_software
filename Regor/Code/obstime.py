import numpy as np
from julian import Night


# wavelength of observation 
lam = 370e-9                        # disk observation (continuous spectrum)
lam1 = 465e-9                       # emission region observation (carbon line spectrum)

# resolution time
delt = 5e-9

# the photon detector quantum efficiency
ita = 0.30

# loss in system
loss = 0.50


# read the starting, ending and number of step of observation
f = open("data/julian.txt","r")

obj = []
N = 3
for i in range(N):
    l = f.readline().split()
    l = [float(s) for s in l]
    obj.append(l)

start = obj[0]                           # list of starting time of observation
end = obj[1]                             # list of ending time of observation
step = obj[2]                            # list of step have been taken for average time 30 minutes
step = [int(a) for a in step]

julday = []
for i in range(0, Night, 1): 
    j = step[i]
    jday = np.linspace(start[i], end[i], step[i])
    julday.append(jday)
               
jd = np.concatenate(julday)             # total observation time according to Julian time

