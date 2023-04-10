import numpy as np
import matplotlib.pyplot as plt
from julian import step
from obspos import rad, x, y
from obstime import lam, lam1
from model import hbt
from starpara import X, Y, t, dx1, dy1, zs
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

# plot the residuals
plt.plot(dx1, dy1, '-', color='darkred')
plt.title('The residuals for 17 day of observation')
plt.savefig(rfig + "resi" + ".png")
plt.show()

# Compare the fitted orbit
plt.plot(X, Y, '.-.', color='indigo', label='simulated orbit')
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
    plt.plot(jds, gts, '.-.', color='green')
    plt.plot(jds, gtsn, '.', label = (east, north))
    plt.errorbar(jds, gtsn, yerr=std(gts), color='blue', fmt='o')
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('correlated Signal $g_w(u)$', fontsize=12)
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.legend()
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
    g = hbt(gX,gY,ts,lam1,lam1)                        # lam=lam1 for 465nm line
    plt.contourf(gX,gY,abs(g))
    plt.plot(xts,yts, '.', color='magenta')
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig(dfig + "cntr/cntr" + str(i+1) + ".png")
    plt.show()



