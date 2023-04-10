import numpy as np
from julian import dt, nstep

# path to read and save file
drk = 'data/orbit/'

# number of objects in system
fil = open(drk + "para.txt")                             # open the txt file with object's parameter
fil.readline()
N = sum(1 for line in fil)                               # number of objects in system                     
fil.close()

# calculate total energy and ratio of Kin and pot energy
mass = np.zeros(N)
pos = np.zeros(shape=(N,3))
vel = np.zeros(shape=(N,3))
def checks():
    kin = pot = 0
    for i in range(N):
        kin += mass[i]*np.sum(vel[i]*vel[i])/2           # kinetic energy for velocity is three dimensional
        for j in range(i):
            dr = pos[i]-pos[j]
            r = np.sum(dr*dr)**.5                        # distance between two object
            pot -= mass[i]*mass[j]/r                     # potential energy, masses are in terms of GMsun
    print('E and vir',kin+pot,-kin/pot)                  # total energy and ratio of kintetic and potential energy (virial theorm)

# open the file to read initial position and velocity
fil = open(drk+"ini.txt","r")
l = fil.readline().split()                     
for i in range(len(l)):                             
    mass[i] = float(l[i])                                # mass

for i in range(N):
    l = fil.readline().split()
    for j in range(len(l)):
        pos[i][j] = float(l[j])                          # position

for i in range(N):
    l = fil.readline().split()
    for j in range(len(l)):
        vel[i][j] = float(l[j])                          # velocity

# calculate the center of mass of position and velocity
cmx = 0*pos[0]
cmv = 0*vel[0]
for i in range(N):
    cmx += mass[i]*pos[i]
    cmv += mass[i]*vel[i]

cmx /= np.sum(mass)                                      # center of mass position
cmv /= np.sum(mass)                                      # center of mass velocity
print('cmx :', cmx)
print('cmv :', cmv)
print('At the end of integration the total enery and ratio of kin and pot for N object')

# Barycentric position and velocity
for i in range(N):
    pos[i] -= cmx                                        # position of Nth body w.r.t. center of mass
    vel[i] -= cmv                                        # velocity of Nth body w.r.t. center of mass velocity

# calculate the next position and velocity (z component) using Leapfrog algorithm
x = np.zeros((N, nstep), dtype = float)
y = np.zeros((N, nstep), dtype = float)
z = np.zeros((N, nstep), dtype = float)
t_stp = np.zeros((nstep,1), dtype = float)
vz = np.zeros((N, nstep), dtype = float)
for t in range(nstep):
    if t % dt == 0:
         checks()
    pos += vel*dt/2                                      # intrim position from an old position (initial position first)
    acc = 0*pos                                          # the accelerations evaluated at the interim position
    t_stp[t] = t_stp[t-1]+ dt                            # the next time step
    for i in range(N):
        for j in range(i):
            dr = pos[i]-pos[j]                           # respective position difference
            r2 = dr[0]**2 + dr[1]**2 + dr[2]**2          # (i,j,k)
            r3 = r2**1.5                                 # respective position 
            acc[j,:] += mass[i]*dr/r3                    # acceleration 
            acc[i,:] -= mass[j]*dr/r3                    
            
    vel += acc*dt                                        # new velocity after steps
    pos += vel*dt/2                                      # new position after steps

    x[:,t] = pos[:,0]                                    # xth position of all N object
    y[:,t] = pos[:,1]                                    # yth position of all N object
    z[:,t] = pos[:,2]                                    # zth position of all N object

    vz[:,t] = vel[:,2]                                   # velocity along the line of sight of all N object

# save the position of N objects in n-files for n-coordinates also and time and z-velocity of N object
np.save(drk + 'x_coord',x)           
np.save(drk + 'y_coord',y)           
np.save(drk + 'z_coord',z)           
np.save(drk + 't_coord',t_stp)  
np.save(drk + 'vz',vz)            

# load all coordinates file of position to create orbit position in x, y and (x,y)
x = np.load(drk + 'x_coord.npy')
y = np.load(drk + 'y_coord.npy')
z = np.load(drk + 'z_coord.npy')
t = np.load(drk + 't_coord.npy')

slice = 1
xs_0 = x[0,::slice]                                      # x-coordinate of Blue-supergiant               
ys_0 = y[0,::slice]                                      # y-coordinate of Blue-supergiant
xs_1 = x[1,::slice]                                      # x-coordinate of Wolf-Rayet           
ys_1 = y[1,::slice]                                      # y-coordinate of Wolf-Rayet                         
ts = t[::slice]

dxAB = dyAB = np.zeros((N, nstep), dtype = float)
dAB = np.zeros((2,len(xs_0)), dtype = float)
dAB[0,:] = xs_0 - xs_1            
dAB[1,:] = ys_0 - ys_1            
dxAB = x[0,::] - x[1,::]         
dyAB = y[0,::] - y[1,::]          
np.save(drk + 'dxAB', dxAB)                               # x-coordinate of binary ow        
np.save(drk + 'dyAB', dyAB)                               # y coordinate of binary ow
np.save(drk + 'dAB', dAB)                                 # (x,y) coordinate of binary ow

