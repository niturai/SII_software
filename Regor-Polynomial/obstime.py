import numpy as np

# path to read file
drk = 'data/intfer/'

# wavelength of observation 
lam = 370e-9                            # disk observation (continuous spectrum)
lam1 = 465e-9                           # emission region observation (carbon line spectrum)

# resolution time
delt = 5e-9

# the photon detector quantum efficiency
ita = 0.30

# loss in system
loss = 0.50

# read the starting, ending and number of step of observation
f = open(drk + "julian.txt","r")
obj = []
N = 3
for i in range(N):
    l = f.readline().split()
    l = [float(s) for s in l]
    obj.append(l)

start = obj[0]                           # list of starting julian time of observation
end = obj[1]                             # list of ending julian time of observation
step = obj[2]                            # list of step have been taken for average time
step = [int(a) for a in step]

# the julian day of observation 
julday = []
for i in range(0, len(step), 1): 
    j = step[i]
    jday = np.linspace(start[i], end[i], step[i])
    julday.append(jday)               
jd = np.concatenate(julday)              # total observation time according to Julian time

