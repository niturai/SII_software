import numpy as np
import matplotlib.pyplot as plt
from obslen import base, lam, lam1
from intfer2 import hbt, X, Y, gX, gY

plt.figure(figsize=(8,8))


# observation Julian day and signals
jd = np.load("data_final/jd.npy")
xt = np.load("data_final/xt.npy") 
yt = np.load("data_final/yt.npy")
gt = np.load("data_final/sig.npy")
gtn = np.load("data_final/sig_noise.npy")
gt1 = np.load("data_final/sig1.npy")
gtn1 = np.load("data_final/sig_noise1.npy")


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
plt.savefig("data_final/fig/sgn_wrl.png")
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
plt.savefig("data_final/fig/sgnoise_wrl.png")
plt.show()    


# plot signal without noise
plt.plot(jd,gt1[0],color='blue', label = (east[0], north[0]))
plt.plot(jd,gt1[1],color='green', label = (east[1], north[1]))
plt.plot(jd,gt1[2],color='red', label = (east[2], north[2]))
plt.plot(jd,gt1[3],color='cyan', label = (east[3], north[3]))
plt.plot(jd,gt1[4],color='purple', label = (east[4], north[4]))
plt.plot(jd,gt1[5],color='olive', label = (east[5], north[5]))
plt.xlabel('Julian Day', ha='center', fontsize=12)
plt.ylabel('correlated Signal $g_w(u)$', fontsize=12)
plt.legend()
plt.savefig("data_final/fig/sgn_rl.png")
plt.show()

# plot signal with noise
plt.plot(jd,gtn1[0],color='blue', label = (east[0], north[0]))
plt.plot(jd,gtn1[1],color='green', label = (east[1], north[1]))
plt.plot(jd,gtn1[2],color='red', label = (east[2], north[2]))
plt.plot(jd,gtn1[3],color='cyan', label = (east[3], north[3]))
plt.plot(jd,gtn1[4],color='purple', label = (east[4], north[4]))
plt.plot(jd,gtn1[5],color='olive', label = (east[5], north[5]))
plt.xlabel('Julian Day', ha='center', fontsize=12)
plt.ylabel('correlated Signal $g_w(u)$ with noise', fontsize=12)
plt.legend()
plt.savefig("data_final/fig/sgnoise_rl.png")
plt.show()    


# compare the signal between roche lobe and without roche lobe
text_kwargs = dict(ha='center', va='center', fontsize=14)
for i in range(len(base)):
    plt.plot(jd,abs(gt[i]),color='green', label = "without roche lobe" )
    plt.plot(jd,abs(gt1[i]),color='red', label = "with roche lobe")
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('correlated Signal $g_w(u)$', fontsize=12)
    plt.title('baseline (E,N) = ' + '(' + str(east[i]) + ',' + str(north[i]) + ')')
    plt.legend()
    plt.savefig("data_final/fig/comp" + str(i) + ".png")
    plt.show()

# compare the noisy signal between roche lobe and without roche lobe
text_kwargs = dict(ha='center', va='center', fontsize=14)
for i in range(len(base)):
    plt.plot(jd,abs(gtn[i]),color='green', label = "without roche lobe")
    plt.plot(jd,abs(gtn1[i]),color='red', label = "with roche lobe")
    plt.xlabel('Julian Day', ha='center', fontsize=12)
    plt.ylabel('correlated Signal $g_w(u)$ with noise', fontsize=12)
    plt.title('baseline (E,N) = ' + '(' + str(east[i]) + ',' + str(north[i]) + ')')
    plt.legend()
    plt.savefig("data_final/fig/comp_noise" + str(i) + ".png")
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
plt.savefig("data_final/fig/track.png")
plt.show()

g = hbt(gX,gY,X,Y,lam,lam1)
g1 = hbt(gX,gY,X,Y,lam1,lam1)


# plot the signal for each baseline with rotation of earth
text_kwargs = dict(ha='center', va='center', fontsize=14)
for i in range(len(base)):
    plt.contour(gX[i],gY[i],abs(g[i]),colors='green', linestyles='dashed')
    plt.contour(gX[i],gY[i],abs(g1[i]),colors='red', linestyles='dotted')
    plt.plot(xt[i],yt[i],'.',color='cyan')
    plt.title('baseline (E,N) = ' + '(' + str(east[i]) + ',' + str(north[i]) + ')')
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.savefig("data_final/fig/cntr_comp" + str(i) + ".png")
    plt.show()


