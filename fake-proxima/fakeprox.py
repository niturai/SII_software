import numpy as np
from scipy.optimize import brentq
from fourier import grids
from sources import blob, smooth
import matplotlib.pyplot as pl


def teff(N,ds,lam):
    ch = pl.imread('eklein.png')
    sx,sy,x,y = grids(ds,N,lam)
    Tmap = 0*sx
    Tmap[N//2-100:N//2+100,N//2-100:N//2+100] = ch[::-1,:,0]
    T = Tmap.reshape(N*N)
    px = T > np.max(T)/10
    Ts = 0*T
    Ts[px] = 3000
    vari = 0*T
    vari[px] = T[px]/np.mean(T[px])
    Tmap = Ts.reshape((N,N))    
    vari = vari.reshape((N,N))
    return sx,sy,x,y,Tmap,vari

