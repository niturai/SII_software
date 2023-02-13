import numpy as np
import matplotlib.pyplot as plt
from obslen import lam, lam1
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
plt.plot(jd,gt, '.-.', color='blue', label = (east, north))
plt.xlabel('Julian Day', ha='center', fontsize=12)
plt.ylabel('correlated Signal $g_w(u)$', fontsize=12)
plt.legend()
plt.savefig("fig/signal.png")
plt.show()

# plot signal with noise
plt.plot(jd,gtn, '.-.', color='blue', label = (east, north))
plt.xlabel('Julian Day', ha='center', fontsize=12)
plt.ylabel('Noisy correlated Signal $g_w(u)$', fontsize=12)
plt.legend()
plt.savefig("fig/sgnoise.png")
plt.show()    


g = hbt(gX,gY,X,Y,lam)

# plot the signal for each baseline with rotation of earth
text_kwargs = dict(ha='center', va='center', fontsize=14)
plt.contourf(gX,gY,abs(g))
plt.plot(xt,yt, '.-.', color='cyan')
plt.title('baseline (E,N) = ' + '(' + str(east) + ',' + str(north) + ')')
plt.colorbar()
plt.gca().set_aspect('equal')
plt.savefig("fig/cntr_comp" + ".png")
plt.show()

