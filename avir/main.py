import numpy as np
import matplotlib.pyplot as pl
from visib import sbright, correldens
from graphics import draw

from spica import teff
from layout import Layout

pl.style.use('dark_background')

source = Layout()
source.set_radec((13,25,11.579),(-11,9,40.75))
source.set_latlon((28,45,25.79),(-17,53,14.39)) # CTA-N + MAGIC

def run():
    F = 40
    N = 1024
    ds = 1e-10
    lam = 400e-9
    for fr in range(F):
        jd = 246e4 + fr/10
        sx,sy,x,y,Tmap = teff(N,ds,lam,jd)
        til = ('Teff on JD = %10.2f' % jd)
        draw(sx,sy,Tmap,8,'sky',cmap='gray',title=til)
        #pl.pause(1)
        pl.savefig('tmp/sky'+'%i'%fr)
        S = sbright(Tmap,lam,ds)
        f = correldens(S,lam)
        #f *= 20 * (1e9 * 300)**.5  # A = 20 m^2 and t_obs = 5 min
        f /= np.max(f)
        til = ('normalized correlation at %i nm' % (1e9*lam+0.5))
        draw(x,y,f,16,'ground',ceil=1,cmap='inferno',title=til)
        u,v,w = source.get_uvw(jd,dx,dy,0*dx)
        hx,hy,hz = source.get_xyz(jd,0,0,1)
        print('time %.3f altitude %.3f' % (jd-246e4,hz))
        if hz > 0:
            pl.plot(u,v,'+',color='white')
        #pl.pause(1)
        pl.savefig('tmp/cdens%i'%fr)
        print(fr)

diam_mst = 12
diam_magic = 17
diam_lst = 23

flags = 'MAGIC'
x = []
y = []
a = []
fil = open('D25.txt')
while True:
    s = fil.readline()
    if not s:
        break
    xs,ys,c = s.split()
    if c == 'MAGIC':
        if 'MAGIC' not in flags:
            continue
        a.append(diam_magic)
    elif c == 'LST':
        if 'LST' not in flags:
            continue
        a.append(diam_lst)
    elif c == 'MST' or c == 'ZST':
        if 'MST' not in flags:
            continue
        a.append(diam_mst)
    print(c)
    x.append(4*float(xs))
    y.append(-4*float(ys))
    


dx = []
dy = []
for i in range(len(x)):
    for j in range(len(x)):
        if i != j:
            dx.append(x[i]-x[j])
            dy.append(y[i]-y[j])
dx = np.array(dx)
dy = np.array(dy)

for i in range(len(dx)):
    print(dx[i],dy[i])

run()

