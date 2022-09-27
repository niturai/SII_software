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


#Orbit plot of gamma 2 velorum
plt.close()
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = [8,8]
plt.plot(x[0,:], y[0,:],color = 'red', linestyle= 'dashdot', linewidth = 1.2, label = 'Gamma  2A')
plt.plot(x[1,:], y[1,:],color = 'green', linestyle= 'dotted', linewidth = 1.2, label = 'Gamma 2B')
plt.title('Orbit plots', pad = 10)
plt.xlabel('x-Coordinates in ls', fontsize=18, color='black')
plt.ylabel('y-coordinate in ls', fontsize=18, color='black')
plt.legend()
fname = 'data/fig/orbit' + '.png'
plt.savefig(fname)
plt.show()


# orbit plot of binary gamma2 velorum
plt.plot(dxAB, dyAB, color = 'red', linestyle= 'dashdot', linewidth = 1.2, label = 'Gamma 2')
plt.title('Relative orbit plots', pad = 10)
plt.xlabel('x-Coordinates in ls', fontsize=18, color='black')
plt.ylabel('y-coordinate in ls', fontsize=18, color='black')
plt.legend()
fname = 'data/fig/Rltvobt' + '.png'
plt.savefig(fname)
plt.show()


# orbit plot of binary gamma 1 velorum
plt.close()
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = [8,8]
plt.plot(x[2,:], y[2,:],color = 'blue', linestyle= 'dotted', linewidth = 1.2, label = 'Gamma  1')
plt.title('Orbit plots', pad = 10)
plt.xlabel('x-Coordinates in ls', fontsize=18, color='black')
plt.ylabel('y-coordinate in ls', fontsize=18, color='black')
plt.legend()
fname = 'data/fig/orbit2' + '.png'
plt.savefig(fname)
plt.show()


# orbit plot of all object
plt.close()
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = [8,8]
plt.plot(x[0,:], y[0,:],color = 'red', linestyle= 'dashed', linewidth = 1.2, label = 'Gamma  2A')
plt.plot(x[1,:], y[1,:],color = 'green', linestyle= 'dotted', linewidth = 1.2, label = 'Gamma  2B')
plt.plot(x[2,:], y[2,:],color = 'blue', linestyle= 'dotted', linewidth = 1.2, label = 'Gamma  1')
plt.title('Orbit plots', pad = 10)
plt.xlabel('x-Coordinates in ls', fontsize=18, color='black')
plt.ylabel('y-coordinate in ls', fontsize=18, color='black')
plt.legend()
fname = 'data/fig/orbit3' + '.png'
plt.savefig(fname)
plt.show()
