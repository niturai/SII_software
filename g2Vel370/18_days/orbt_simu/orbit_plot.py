import numpy as np
import matplotlib.pyplot as plt

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


# compare orbit with 18 days
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = [10,10]
plt.plot(-x[0,:], -y[0,:],color = 'red', linestyle= 'dashdot', linewidth = 1.2, label = 'Blue Supergiant')
plt.plot(-x[1,:], -y[1,:],color = 'green', linestyle= 'dotted', linewidth = 1.2, label = 'Wolf Rayet')
plt.plot(-x[0,:865], -y[0,:865],color = 'blue', linewidth = 1.2, label = 'BS for 18 days')
plt.plot(-x[1,:865], -y[1,:865],color = 'blue', linewidth = 1.2, label = 'WR for 18 days')
plt.plot(-x[0,0],-y[0,0],marker='*', color='black',markersize=6, label = 'BS initial point' )
plt.plot(-x[1,0],-y[1,0],marker='*', color='black',markersize=6, label = 'WR initial point')
plt.xlabel('X in light second')
plt.ylabel('Y in light second')
plt.title("Orbit of Blue Supergiant and Wolf Rayet")
plt.legend()
plt.savefig('fig/pos18.png')
plt.show()


plt.plot(-dAB[0,:],-dAB[1,:], color='cyan',linestyle= 'dashdot', linewidth = 1.2, label = 'orbit of Regor')
plt.plot(-dAB[0,:865],-dAB[1,:865], color='blue', linewidth = 1.2, label = 'Regor for 18 days')
plt.plot(-dAB[0,0],-dAB[1,0],marker='*',color='black', markersize=6, label='initial point')
plt.xlabel('X in light second')
plt.ylabel('Y in light second')
plt.title("Orbit of Regor")
plt.legend()
plt.savefig('fig/orb18.png')
plt.show()

