import numpy as np
import matplotlib.pyplot as plt
from julian import step
from obspos import rad, x, y
from obstime import lam, lam1
from model import hbt
from starpara import X, Y
from rotation import grids

# observation Julian day and signals
jd = np.load("data/370/jd.npy")
xt = np.load("data/370/xt.npy") 
yt = np.load("data/370/yt.npy")
gt = np.load("data/370/sig.npy")
gtn = np.load("data/370/sig_noise.npy")
gt1 = np.load("data/465/sig.npy")
gtn1 = np.load("data/465/sig_noise.npy")

# calculate the standard deviation
def std(gt):
    return np.std(gt)

# three dimensional baseline 
east = x
north = y

# calculate the steps from data of observation for 17 days
st = np.cumsum(step)

# plot the signal for 370nm
for i in range(0, len(st)):
    if i==0:
       jds = jd[ : st[i]]
       gts = gt[ : st[i]]
       gtsn = gtn[ : st[i]]
    else:
       jds = jd[st[i-1]:st[i]]
       gts = gt[st[i-1]:st[i]]
       gtsn = gtn[st[i-1]:st[i]]
    plt.plot(jds, gts, '.-.', color='green')
    plt.plot(jds, gtsn, '.', label = (east, north))
    plt.errorbar(jds, gtsn, yerr=std(gts), color='blue', fmt='o')
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('correlated Signal $g_w(u)$ for 370nm', fontsize=12)
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig("fig/370/sig/sig" + str(i+1) + ".png")
    plt.show()


# plot the signal for 465nm
for i in range(0, len(step)):
    if i==0:
       jds = jd[ : st[i]]
       gts = gt1[ : st[i]]
       gtsn = gtn1[ : st[i]]
    else:
       jds = jd[st[i-1]:st[i]]
       gts = gt1[st[i-1]:st[i]]
       gtsn = gtn1[st[i-1]:st[i]]
    plt.plot(jds, gts, '.-.', color='green')
    plt.plot(jds, gtsn, '.', label = (east, north))
    plt.errorbar(jds, gtsn, yerr=std(gts), color='blue', fmt='o')
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('correlated Signal $g_w(u)$ for 465nm', fontsize=12)
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig("fig/465/sig/sig" + str(i+1) + ".png")
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
    g = hbt(gX,gY,Xs,Ys,lam,lam1)
    plt.contourf(gX,gY,abs(g))
    plt.plot(xts,yts, '.', color='magenta')
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig("fig/370/cntr/cntr" + str(i+1) + ".png")
    plt.show()

# Signal for 465nm according to roation of earth
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
    g = hbt(gX,gY,Xs,Ys,lam1,lam1)
    plt.contourf(gX,gY,abs(g))
    plt.plot(xts,yts, '.', color='magenta')
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig("fig/465/cntr/cntr" + str(i+1) + ".png")
    plt.show()

