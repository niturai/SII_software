import numpy as np
from scipy.optimize import brentq
from fourier import grids
from sources import blob, smooth

pi, cos, sin = np.pi, np.cos, np.sin

def pos(tph,e,omega,I,Omega):
    ec = (1-e*e)**.5
    psi = brentq(lambda psi: psi - e*sin(psi) - tph, 0, 2*pi)
    cs,sn = cos(psi), sin(psi)
    x,y = cs-e, ec*sn
    cs,sn = cos(omega), sin(omega)
    x,y = x*cs - y*sn, x*sn + y*cs
    x,y = x, y*cos(I)
    cs,sn = cos(Omega), sin(Omega)
    x,y = x*cs - y*sn, x*sn + y*cs
    return x,y


def teff(N,ds,lam,jd):
    Rsun = 2.3
    pc = 1.03e8
    dis = 80*pc
    R1, T1 = 6.4*Rsun/dis, 25e3
    R2, T2 = 3.2*Rsun/dis, 20e3
    a = 28*Rsun/dis
    q = 0.63
    a1 = q/(1+q)*a
    a2 = -1/(1+q)*a
    e = 0.14
    omega = 140*pi/180
    I = (180-63)*pi/180
    Omega = 130*pi/180
    jdperi = 2440678.09
    Pinv = 0.249091
    tph = 2*pi*(jd-jdperi)*Pinv
    tph = tph % (2*pi)
    xs,ys = pos(tph,e,omega,I,Omega)
    sx,sy,x,y = grids(ds,N,lam)
    Tmap1 = T1 * blob(sx,sy,R1,a1*xs,a1*ys)
    Tmap2 = T2 * blob(sx,sy,R2,a2*xs,a2*ys)
    Tmap = smooth(Tmap1+Tmap2,3)
    return sx,sy,x,y,Tmap

