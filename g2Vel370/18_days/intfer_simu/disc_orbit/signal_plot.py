import numpy as np
import matplotlib.pyplot as plt
from obslen import Night, lam, lam1
from intfer import x, y, hbt, X, Y, gX, gY


# observation Julian day and signals
jd = np.load("data/jd.npy")
xt = np.load("data/xt.npy") 
yt = np.load("data/yt.npy")
gt = np.load("data/sig.npy")
gtn = np.load("data/sig_noise.npy")


# three dimensional baseline 
east = x
north = y

# plot signal without noise
for i in range(0, Night, 1):
    plt.plot(jd[i],gt[i], '.-.', color='blue', label = (east, north))
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('correlated Signal $g_w(u)$', fontsize=12)
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig("fig/signal" + str(i+1) + ".png")
    plt.show()

# plot signal with noise
for i in range(0, Night, 1):
    plt.plot(jd[i],gtn[i], '.-.', color='blue', label = (east, north))
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('Noisy correlated Signal $g_w(u)$', fontsize=12)
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig("fig/sig_noise" + str(i+1) + ".png")
    plt.show()    



# plot the signal for each baseline with rotation of earth
for i in range(0, Night, 1):
    g = hbt(gX[i,:],gY[i,:],X,Y,lam)
    plt.contourf(gX[i],gY[i],abs(g))
    plt.plot(xt[i],yt[i], '.', color='black')
    plt.title(str(i+1) + 'th' + ' day of observation')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig("fig/cntr" + str(i+1) + ".png")
    plt.show()


