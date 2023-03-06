import numpy as np
import matplotlib.pyplot as plt
from julian import Night, dst, step, den

#********************************************#
# load dAB file
dAB = np.load("data/dAB.npy")               # (x,y) of gamma2 velorum
dxAB = np.load("data/dxAB.npy")             # x of gamma 2 velorum's star
dyAB = np.load("data/dyAB.npy")             # y of gamma 2 velorum's srar

# load coordinates file of all objects
x = np.load("data/x_coord.npy")
y = np.load("data/y_coord.npy")
z = np.load("data/z_coord.npy")
t = np.load("data/t_coord.npy")

# compare orbit of Blue Supergiant and Wolf Rayet with 17 days
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = [10,10]
plt.plot(-x[0,:], -y[0,:],color = 'red', linestyle= 'dashdot', linewidth = 3)
plt.plot(-x[1,:], -y[1,:],color = 'green', linestyle= 'dotted', linewidth = 3)
plt.plot(-x[0,:den[16]], -y[0,:den[16]],color = 'black', linewidth = 3, label = 'BS for 17 days')
plt.plot(-x[1,:den[16]], -y[1,:den[16]],color = 'blue', linewidth = 3, label = 'WR for 17 days')
plt.plot(-x[0,den[16]-1:den[16]], -y[0,den[16]-1:den[16]],marker='o',color = 'indigo',markersize=6)
plt.plot(-x[1,den[16]-1:den[16]], -y[1,den[16]-1:den[16]],marker='o',color = 'purple',markersize=6)
plt.plot(-x[0,0],-y[0,0],marker='x', color='black',markersize=10, label = 'initial point' )
plt.plot(-x[1,0],-y[1,0],marker='x', color='black',markersize=10)
plt.xlabel('X in light second')
plt.ylabel('Y in light second')
plt.title("Orbit of Blue Supergiant and Wolf Rayet")
plt.legend()
plt.savefig('fig/orbit/pos17.png')
plt.show()

# compare orbit of Regor with 16 days
plt.plot(-dAB[0,:],-dAB[1,:], color='cyan',linestyle= 'dashdot', linewidth = 3)
plt.plot(-dAB[0,:den[16]],-dAB[1,:den[16]], color='blue', linewidth = 3, label = 'Regor for 17 days')
plt.plot(-dAB[0,den[16]-1:den[16]],-dAB[1,den[16]-1:den[16]], marker='o', color='purple', markersize=6)
plt.plot(-dAB[0,0],-dAB[1,0],marker='x',color='black', markersize=10, label='initial point')
plt.xlabel('X in light second')
plt.ylabel('Y in light second')
plt.title("Orbit of Regor")
plt.legend()
plt.savefig('fig/orbit/orb17.png')
plt.show()


# the position of Blue Supergiant and Wolf Rayet, on the observation day
for i in range(0,Night,1):
    plt.plot(-x[0, :], -y[0, :],color = 'red', linestyle= 'dashdot', linewidth = 3)
    plt.plot(-x[1, :], -y[1, :],color = 'green', linestyle= 'dotted', linewidth = 3)
    plt.plot(-x[0, dst[i] : den[i]], -y[0, dst[i] : den[i]],color = 'cyan', linewidth = 3, label = 'BS on the observation time')
    plt.plot(-x[1, dst[i] : den[i]], -y[1, dst[i] : den[i]],color = 'blue', linewidth = 3, label = 'WR on the observation time')
    plt.plot(-x[0,0],-y[0,0],marker='x', color='black',markersize=10, label = 'initial point' )
    plt.plot(-x[1,0],-y[1,0],marker='x', color='black',markersize=10)
    plt.xlabel('X in light second')
    plt.ylabel('Y in light second')
    plt.title("Blue Supergiant and Wolf Rayet " + str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig('fig/orbit/pos/pos' + str(i+1) + '.png')
    plt.show()

# the orbit of Regor on observation day
for i in range(0,Night,1):
    j = step[i]
    plt.plot(-dAB[0, :],-dAB[1, :], color='cyan',linestyle= 'dashdot', linewidth = 3)
    plt.plot(-dAB[0, dst[i] : den[i]],-dAB[1, dst[i] : den[i]], color='blue', linewidth = 3, label = 'Regor on the observation time')
    plt.plot(-dAB[0,0],-dAB[1,0],marker='x',color='black', markersize=10, label='initial point')
    plt.xlabel('X in light second')
    plt.ylabel('Y in light second')
    plt.title("Regor " + str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig('fig/orbit/orb/orb' + str(i+1) + '.png')
    plt.show()

