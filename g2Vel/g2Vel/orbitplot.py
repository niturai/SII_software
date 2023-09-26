import numpy as np
import matplotlib.pyplot as plt
from julian import Night, dst, step, den

# path to read the file
drk = 'data/orbit/'

# path to save the file
dfig = 'fig/orbit/'

#********************************************#
# load dAB file
dAB = np.load(drk + "dAB.npy")               # (x,y) of gamma2 velorum
dxAB = np.load(drk + "dxAB.npy")             # x of gamma 2 velorum's star
dyAB = np.load(drk + "dyAB.npy")             # y of gamma 2 velorum's srar

# load coordinates file of all objects
x = np.load(drk + "x_coord.npy")
y = np.load(drk + "y_coord.npy")
z = np.load(drk + "z_coord.npy")
t = np.load(drk + "t_coord.npy")
vz = np.load(drk + "vz.npy")

# plot all the simulation
plt.rcParams.update({'font.size': 10})
plt.rcParams["figure.figsize"] = [8,8]

# the velocity of Blue Supergiant and Wolf Rayet
plt.plot(t,vz[0,:], color='red', label='Blue Supergiant')
plt.plot(t,vz[1,:], color='green', label='Wolf-Rayet')
plt.xlabel('time in seconds')
plt.ylabel('Velovity in ls per sec')
plt.title('Compare the velocities of Regor objects')
plt.legend()
plt.savefig(dfig + 'vel.png')
plt.show()

# compare orbit of Blue Supergiant and Wolf Rayet with 17 days
plt.plot(-x[0,:], -y[0,:],color = 'red', linestyle= 'dashdot', linewidth = 3)
plt.plot(-x[1,:], -y[1,:],color = 'green', linestyle= 'dotted', linewidth = 3)
plt.plot(-x[0,:den[16]], -y[0,:den[16]],color = 'black', linewidth = 3, label = 'BS for 17 days')
plt.plot(-x[1,:den[16]], -y[1,:den[16]],color = 'blue', linewidth = 3, label = 'WR for 17 days')
plt.plot(-x[0,den[16]:den[16]+1], -y[0,den[16]:den[16]+1],marker='o',color = 'indigo',markersize=6)
plt.plot(-x[1,den[16]:den[16]+1], -y[1,den[16]:den[16]+1],marker='o',color = 'purple',markersize=6)
plt.plot(-x[0,0],-y[0,0],marker='x', color='black',markersize=8, label = 'BS initial point' )
plt.plot(-x[1,0],-y[1,0],marker='x', color='blue',markersize=8, label = 'WR initial point')
plt.xlabel('X in light second')
plt.ylabel('Y in light second')
plt.title("Orbit of Blue Supergiant and Wolf Rayet")
plt.legend()
plt.savefig(dfig + 'pos17.png')
plt.show()

# compare orbit of Regor with 16 days
plt.plot(-dAB[0,:],-dAB[1,:], color='cyan',linestyle= 'dashdot', linewidth = 3)
plt.plot(-dAB[0,:den[16]],-dAB[1,:den[16]], color='fuchsia', linewidth = 3, label = 'Regor for 17 days')
plt.plot(-dAB[0,den[16]:den[16]+1],-dAB[1,den[16]:den[16]+1], marker='o', color='blue', markersize=6)
plt.plot(-dAB[0,0],-dAB[1,0],marker='x',color='black', markersize=8, label='initial point')
plt.xlabel('X in light second')
plt.ylabel('Y in light second')
plt.title("Orbit of Regor")
plt.legend()
plt.savefig(dfig + 'orb17.png')
plt.show()


# the position of Blue Supergiant and Wolf Rayet, on the observation day
for i in range(0,Night,1):
    plt.plot(-x[0, :], -y[0, :],color = 'red', linestyle= 'dashdot', linewidth = 3)
    plt.plot(-x[1, :], -y[1, :],color = 'green', linestyle= 'dotted', linewidth = 3)
    plt.plot(-x[0, dst[i] : den[i]], -y[0, dst[i] : den[i]],color = 'black', linewidth = 3, label = 'BS on the observation time')
    plt.plot(-x[1, dst[i] : den[i]], -y[1, dst[i] : den[i]],color = 'blue', linewidth = 3, label = 'WR on the observation time')
    plt.plot(-x[0,0],-y[0,0],marker='o', color='red',markersize=8, label = 'BS initial point' )
    plt.plot(-x[1,0],-y[1,0],marker='o', color='green',markersize=8, label = 'WR initial point')
    plt.xlabel('X in light second')
    plt.ylabel('Y in light second')
    plt.title("Blue Supergiant and Wolf Rayet " + str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig(dfig + 'pos/pos' + str(i+1) + '.png')
    plt.show()

# the orbit of Regor on observation day
for i in range(0,Night,1):
    j = step[i]
    plt.plot(-dAB[0, :],-dAB[1, :], color='cyan',linestyle= 'dashdot', linewidth = 3)
    plt.plot(-dAB[0, dst[i] : den[i]],-dAB[1, dst[i] : den[i]], color='fuchsia', linewidth = 3, label = 'Regor on the observation time')
    plt.plot(-dAB[0,0],-dAB[1,0],marker='o',color='black', markersize=8, label='initial point')
    plt.xlabel('X in light second')
    plt.ylabel('Y in light second')
    plt.title("Regor " + str(i+1) + 'th' + ' day of observation')
    plt.legend()
    plt.savefig(dfig + 'orb/orb' + str(i+1) + '.png')
    plt.show()

