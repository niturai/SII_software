import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from julian import step
from obspos import rad, x, y
from obstime import lam, lam1
from model import hbt
from starpara import R_a, R_b, R1, X, Y
from rotation import grids, gX, gY

# path to read file
drk = 'data/'

# path to save file 
dfig = 'fig/'        # the signal

# observation Julian day and signals
jd = np.load(drk + "jd.npy")
xt = np.load(drk + "xt.npy") 
yt = np.load(drk + "yt.npy")
gt = np.load(drk + "sig.npy")
gtn = np.load(drk + "sig_noise.npy")

# calculate the standard deviation
def std(gt):
    return np.std(gt)

# three dimensional baseline 
east = x
north = y

# calculate the steps from data of observation for 17 days
st = np.cumsum(step)

# plot the signal for given wavelengths
for i in range(0, len(st)):
    if i==0:
       jds = jd[ : st[i]]
       gts = gt[ : st[i]]
       gtsn = gtn[ : st[i]]
    else:
       jds = jd[st[i-1]:st[i]]
       gts = gt[st[i-1]:st[i]]
       gtsn = gtn[st[i-1]:st[i]]
    #plt.plot(jds, gts, '.-.', color='green')
    plt.plot(jds, gtsn, '.', label = (east, north))
    plt.errorbar(jds, gtsn, yerr=std(gts), color='green', fmt='o')
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('correlated Signal $g_w(u)$ with noise', fontsize=12)
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.savefig(dfig + "sig/sig" + str(i+1) + ".png")
    plt.show()


# Signal for 370nm according to roation of earth
for i in range(0, len(step)):
    if i==0:
       Xs = X[ : st[i]]
       Ys = Y[ : st[i]]
       xts = xt[ : st[i]]
       yts = yt[ : st[i]]
    else:
       Xs = X[st[i-1]:st[i]]
       Ys = Y[st[i-1]:st[i]]
       xts = xt[st[i-1]:st[i]]
       yts = yt[st[i-1]:st[i]]
    gX, gY = grids(xts,yts,rad,step[i])
    g = hbt(gX, gY, Xs, Ys, lam1, lam1, R_a, R_b, R1)         # lam=lam1 for 465nm line and R1 for roche lobe
    plt.contourf(gX, gY, abs(g), cmap="Greens")
    plt.plot(xts, yts, '.', color='red')
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig(dfig + "cntr/cntr" + str(i+1) + ".png")
    plt.show()

