import numpy as np
import matplotlib.pyplot as plt
from obslen import base, lam
from intfer import hbt, X, Y, gX, gY

plt.figure(figsize=(8,8))


# observation Julian day and signals
jd = np.load("data/jd.npy")
xt = np.load("data/xt.npy") 
yt = np.load("data/yt.npy")
gt = np.load("data/sig.npy")
gtn = np.load("data/sig_noise.npy")


# three dimensional baseline 
east = np.real(base)
north = np.imag(base)

# plot signal without noise
plt.plot(jd,gt[0],color='blue', label = (east[0], north[0]))
plt.plot(jd,gt[1],color='green', label = (east[1], north[1]))
plt.plot(jd,gt[2],color='red', label = (east[2], north[2]))
plt.plot(jd,gt[3],color='cyan', label = (east[3], north[3]))
plt.plot(jd,gt[4],color='purple', label = (east[4], north[4]))
plt.plot(jd,gt[5],color='olive', label = (east[5], north[5]))
plt.xlabel('Julian Day', ha='center', fontsize=12)
plt.ylabel('correlated Signal $g_w(u)$', fontsize=12)
plt.legend()
plt.savefig("data/fig/sgn1.png")
plt.show()

# plot signal with noise
plt.plot(jd,gtn[0],color='blue', label = (east[0], north[0]))
plt.plot(jd,gtn[1],color='green', label = (east[1], north[1]))
plt.plot(jd,gtn[2],color='red', label = (east[2], north[2]))
plt.plot(jd,gtn[3],color='cyan', label = (east[3], north[3]))
plt.plot(jd,gtn[4],color='purple', label = (east[4], north[4]))
plt.plot(jd,gtn[5],color='olive', label = (east[5], north[5]))
plt.xlabel('Julian Day', ha='center', fontsize=12)
plt.ylabel('correlated Signal $g_w(u)$ with noise', fontsize=12)
plt.legend()
plt.savefig("data/fig/sgn2.png")
plt.show()    


# to get the track of all baseline
text_kwargs = dict(ha='center', va='center', fontsize=14)
plt.plot(xt[0],yt[0],'.',color='blue')
plt.plot(xt[1],yt[1],'.',color='black')
plt.plot(xt[2],yt[2],'.',color='red')
plt.plot(xt[3],yt[3],'.',color='cyan')
plt.plot(xt[4],yt[4],'.',color='purple')
plt.plot(xt[5],yt[5],'.',color='yellow')
plt.gca().set_aspect('equal')
plt.savefig("data/fig/track.png")
plt.show()


# plot the signal for each baseline with rotation of earth
g = hbt(gX,gY,X,Y,lam)
text_kwargs = dict(ha='center', va='center', fontsize=14)
for i in range(len(g)):
    plt.contourf(gX[i],gY[i],abs(g[i]))
    plt.plot(xt[i],yt[i],'.',color='red')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig("data/fig/track" + str(i) + ".png")
    plt.show()

