import numpy as np
import matplotlib.pyplot as pl
from visib import sbright, correldens
from graphics import draw

from regor import teff
from layout import Layout

pl.style.use('dark_background')

source = Layout()
source.set_radec((8,9,31.95013),(-47,20,11.7108))
# source.set_latlon((31,40,30),(-110,57,7))  # VERITAS
# source.set_latlon((28,45,25.79),(-17,53,14.39)) # CTA-N + MAGIC
# source.set_latlon((32,48,0),(79,0,0)) # Hanle
source.set_latlon((-23,16,25.4856),(16,31,10.2432))  # HESS


def run():
    F = 260
    N = 1024
    ds = .9e-10
    lam = 465e-9
    for fr in range(1,F):
        jd = 246e5 + (fr-1)/24
        phase = 2*np.pi * ((jd/4) % 1.0)
        sx,sy,x,y,Tmap = teff(N,ds,lam,phase)
        til = ('effective temperature')
        draw(sx,sy,Tmap,4,'sky',cmap='gray',title=til)
        pl.pause(1)
        #pl.savefig('tmp/sky'+'%i'%fr)
        S = sbright(Tmap,lam,ds)
        f = correldens(S,lam)
        f *= 1 * (1e9 * 300)**.5  # A = 1 m^2 and t_obs = 5 min
        f /= np.max(f)
        til = ('normalized correlation at %i nm' % (1e9*lam+0.5))
        draw(x,y,f,8,'ground',ceil=1,cmap='inferno',title=til)
        u,v,w = source.get_uvw(jd,dx,dy,0*dx)
        hx,hy,hz = source.get_xyz(jd,0,0,1)
        print('time %.3f altitude %.3f' % (jd-246e5,hz))
        if hz > 0:
            pl.plot(u,v,'+',color='white')
        pl.pause(1)
        #pl.savefig('tmp/cdens%i'%fr)
        print(fr)

x = []
y = []
fil = open('locs-hess.txt')
while True:
    s = fil.readline()
    if not s:
        break
    xs,ys,c = s.split()
    x.append(float(xs))
    y.append(float(ys))

dx = []
dy = []
for i in range(len(x)):
    for j in range(len(x)):
        if i != j:
            dx.append(x[i]-x[j])
            dy.append(y[i]-y[j])
dx = np.array(dx)
dy = np.array(dy)



run()

