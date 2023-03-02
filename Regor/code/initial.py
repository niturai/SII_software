import numpy as np
import scipy.optimize as sciopt
from julian import start

pi = np.pi
cos = np.cos
sin = np.sin


# rotation of orbit of star system
def rotate(x,y,omega,I,Omega):                 
    cs, sn = cos(omega), sin(omega)
    x, y = x*cs - y*sn, x*sn + y*cs
    cs, sn = cos(I), sin(I)
    x,y,z = x, y*cs, y*sn
    Omega += pi/2                         # rotate (north,east)
    cs,sn = cos(Omega), sin(Omega)
    x,y,z = x*cs - y*sn, x*sn + y*cs, z
    return np.array((x,y,z))


# define the position w.r.t. observer
def posvel(P,a,e,omega,I,Omega,mean_an):
    ec = (1-e*e)**.5
    mean_an = mean_an % (2*pi)                                              # mean anomaly 
    psi = sciopt.brentq(lambda psi: psi - e*sin(psi) - mean_an, 0, 2*pi)    # solve Kepler's equation for eccentric anomaly psi in [0, 2*pi]
    cs,sn = cos(psi), sin(psi)
    x,y = cs-e, ec*sn
    vx,vy = -sn/(1-e*cs), ec*cs/(1-e*cs)
    pos = a * rotate(x,y,omega,I,Omega)
    vel = 2*pi * a/P * rotate(vx,vy,omega,I,Omega)
    return pos,vel

def bary(m,x):
    N = len(m)
    cm = 0*m
    cm[0] = m[0]                                           # m[0], the first object mass
    cx = 0*x
    for n in range(1,N):
        cm[n] = cm[n-1] + m[n]                             # sum of mass of N bodies (0,1,2,3....)
    cx = 0*x
    for n in range(1,N):
        x[n] += cx[n-1]                                    # position of nth body w.r.t. center of mass (m[0])
        cx[n] = (cm[n-1]*cx[n-1] + m[n]*x[n])/cm[n];       # postion of center of mass of N bodies
    b = cx[N-1]                                            # position of barycenter of N bodies
    x[:] -= b                                              # array of 3-dimensional position of N bodies w.r.t. to barycenter
    return x


def getposvel(fname,now=float(start[0])):                  # now is the first day and time of observation, April 26th, 17:35 
    fil = open(fname)
    fil.readline()
    star = []
    N = 3
    for n in range(N):
        pars = fil.readline().split()
        pars = [float(s) for s in pars]
        star.append(pars)                                   # list of parameters of N bodies

    c = 299792458
    day = 86400                                             # one day to seconds
    au = 149597870700 / c                                   # AU to seconds (meter/(meter sec^-1))
    deg = pi/180                                            # degree to radian
    GMsun = 4.9254909e-6                                    # mass of Sun (in second) i.e. GMsun/c**3
   
    mass = np.zeros(N)
    pos = np.zeros((N,3))                                   # three dimensional position
    vel = np.zeros((N,3))                                   # three dimensional velocity

    for n in range(N):
        mass[n] = float(star[n][0]) * GMsun                 # mass of Nth body in second
    for n in range(1,N): 
        P,a,e,I,Omega,ep,omega = star[n][1:]                # all parameters of Nth body
        mean_an = (now-ep)/P * 2*pi                         # mean anomaly (in radian) if (now-ep) and p in days
        P *= day                                            # period in second
        a *= au                                             # semi-major axis in second
        I *= deg                                            # inclination in terms of radian
        Omega *= deg                                        # longitude of the ascending node in terms of radian
        omega *= deg                                        # argument of periapsis in terms of radian
        pos[n],vel[n] = posvel(P,a,e,omega,I,Omega,mean_an) # position and velocity
    return mass, bary(mass,pos), bary(mass,vel)             # return mass , position and velocity of Nth body (initial value) 

mass, pos, vel = getposvel('data/para.txt')                  # read the input parameters of N bodies

fil = open('data/ini.txt', 'w')                              # output as a txt file
         
N_body = len(mass)
for l in range(N_body):
    fil.write("%22.15e" %mass[l])                                 # write the masses
    fil.write('\t')                       
fil.write("\n")                         

for l in range(N_body):
    for k in range(3):                                            # three dimensional position
        fil.write("%22.15e" %pos[l,k])        
        fil.write('\t')
    fil.write("\n")                           

for l in range(N_body):
    for k in range(3):                                            # three dimensional velocity
        fil.write("%22.15e" %vel[l,k])        
        fil.write('\t')
    fil.write("\n")                           

fil.close()                                  

print("mass = ",mass)                         
print("pos = ",pos)                           
print("vel = ",vel)                           
print('center of mass position = ',np.sum(mass*pos[:,0])/np.sum(mass))      

