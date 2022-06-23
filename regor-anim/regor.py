import numpy as np
from scipy.optimize import brentq
from fourier import grids
from sources import blob, smooth

pi, cos, sin = np.pi, np.cos, np.sin

def pos(t,e,omega,I,Omega):
    ec = (1-e*e)**.5
    psi = brentq(lambda psi: psi - e*sin(psi) - t, 0, 2*pi)
    cs,sn = cos(psi), sin(psi)
    x,y = cs-e, ec*sn
    cs,sn = cos(omega), sin(omega)
    x,y = x*cs - y*sn, x*sn + y*cs
    x,y = x, y*cos(I)
    cs,sn = cos(Omega), sin(Omega)
    x,y = x*cs - y*sn, x*sn + y*cs
    return x,y


def teff(N,ds,lam,t=0):
    Rsun = 2.3
    pc = 1.03e8
    dis = 380*pc
    R1, T1 = 17*Rsun/dis, 35e3
    R2, T2 = 6*Rsun/dis, 57e3
    a = 600/dis
    q = 9/28
    a1 = q/(1+q)*a
    a2 = -1/(1+q)*a
    e = 0.326
    omega = 140*pi/180
    I = (180-63)*pi/180
    Omega = 130*pi/180
    xs,ys = pos(t,e,omega,I,Omega)
    sx,sy,x,y = grids(ds,N,lam)
    Tmap1 = T1 * blob(sx,sy,R1,a1*xs,a1*ys)
    Tmap2 = T2 * blob(sx,sy,R2,a2*xs,a2*ys)
    Tmap = smooth(Tmap1+Tmap2,3)
    return sx,sy,x,y,Tmap

