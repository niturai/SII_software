import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from julian import step
from obspos import rad, x, y
from obstime import lam, lam1
from model import hbt
from starpara import X, Y, R_a, R_b, R1, t, dx1, dy1, zs
from rotation import grids, gX, gY

# path to read file
drk = 'data/intfer/465/'

# path to save file 
rfig = 'fig/orbit/'             # the residuals and fitted curve
dfig = 'fig/intfer/465/'        # the signal

# plot the fitting line with residuals
def func(zs, t):
    x_t = y_t = 0
    tp = 1
    for p in range(0, len(zs), 2):
        x_t += zs[p] * tp
        y_t += zs[p+1] * tp
        tp *= t
    return x_t, y_t
Xt, Yt = func(zs, t)

# orbit fitting coefficients for X and Y
x0, x1, x2, x3 = zs[0], zs[2], zs[4], zs[6]
y0, y1, y2, y3 = zs[1], zs[3], zs[5], zs[7]

# plot the residuals
plt.plot(dx1, dy1, '-', color='green')
plt.title('The residuals for 17 day of observation')
plt.savefig(rfig + "resi" + ".png")
plt.show()

# Compare the fitted orbit
plt.plot(X, Y, '.-.', color='red', label='simulated orbit')
plt.plot(Xt, Yt, '.', color='green', label='fited orbit')
plt.title('Compare the fiting curve for 17 day of observation')
plt.legend()
plt.savefig(rfig + "fitorb" + ".png")
plt.show()

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
       ts = t[ : st[i]]
       xts = xt[ : st[i]]
       yts = yt[ : st[i]]
    else:
       ts = t[st[i-1]:st[i]]
       xts = xt[st[i-1]:st[i]]
       yts = yt[st[i-1]:st[i]]
    gX, gY = grids(xts,yts,rad,step[i])
    g = hbt(gX, gY, ts, lam1, lam1, R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3)         # lam=lam1 for 465nm line and R1 for roche lobe
    plt.contourf(gX, gY, abs(g), cmap="Greens")
    plt.plot(xts, yts, '.', color='red')
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig(dfig + "cntr/cntr" + str(i+1) + ".png")
    plt.show()

