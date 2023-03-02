import numpy as np
from julian import dt, nstep

N_body = 3
mass = np.zeros(N_body)
pos = np.zeros(shape=(N_body,3))
vel = np.zeros(shape=(N_body,3))

def checks():
    cmx = 0*pos[0]                                       # center of mass of position
    cmv = 0*pos[0]                                       # center of mass of velocity
    for i in range(N_body):
        cmx += mass[i]*pos[i]
        cmv += mass[i]*vel[i]
    print('cmx: ',cmx)
    print('cmv: ',cmv)
    kin = pot = 0
    for i in range(N_body):
        kin += mass[i]*np.sum(vel[i]*vel[i])/2           # kinetic energy
        for j in range(i):
            dr = pos[i]-pos[j]
            r = np.sum(dr*dr)**.5
            pot -= mass[i]*mass[j]/r                    # potential energy
    print('E and vir',kin+pot,-kin/pot)       

checks()

fil = open("data/ini.txt","r")
l = fil.readline().split()                     
 
for i in range(len(l)):                           
    mass[i] = float(l[i])
print('mass = ', mass)

for i in range(N_body):
    l = fil.readline().split()
    for j in range(len(l)):
        pos[i][j] = float(l[j])
print('pos = ', pos)

for i in range(N_body):
    l = fil.readline().split()
    for j in range(len(l)):
        vel[i][j] = float(l[j])
print('vel = ', vel)

cmx = 0*pos[0]
cmv = 0*vel[0]
for i in range(N_body):
    cmx += mass[i]*pos[i]
    cmv += mass[i]*vel[i]

cmx /= np.sum(mass)                   # center of mass position
cmv /= np.sum(mass)                   # center of mass velocity

for i in range(N_body):
    pos[i] -= cmx                     # position of Nth body w.r.t. center of mass
    vel[i] -= cmv                     # velocity of Nth body w.r.t. center of mass velocity

print('cmx :', cmx)
print('cmv :', cmv)

GMsun = 4.9254909e-6

x = np.zeros((N_body, nstep), dtype = float)
y = np.zeros((N_body, nstep), dtype = float)
z = np.zeros((N_body, nstep), dtype = float)
t_stp = np.zeros((nstep,1), dtype = float)

for t in range(nstep):
    if t % dt == 0:
         checks()
    pos += vel*dt/2                       # intrim position from an old position (initial position first)
    acc = 0*pos                           
    t_stp[t] = t_stp[t-1]+ dt             # the next time step
    for i in range(N_body):
        for j in range(i):
            dr = pos[i]-pos[j]                     # respective position difference
            r2 = dr[0]**2 + dr[1]**2 + dr[2]**2    # (i,j,k)
            r3 = r2**1.5                           # respective position 
            acc[j,:] += mass[i]*dr/r3              # 
            acc[i,:] -= mass[j]*dr/r3              # 
            
    vel += acc*dt                                  # new velocity after steps
    pos += vel*dt/2                                # new position after steps

    x[:,t] = pos[:,0]                              # xth position of all N object
    y[:,t] = pos[:,1]                              # yth position of all N object
    z[:,t] = pos[:,2]                              # zth position of all N object

# save all coordinates in file
np.save('data/x_coord',x)           
np.save('data/y_coord',y)           
np.save('data/z_coord',z)           
np.save('data/t_coord',t_stp)              

# load all coordinates file 
x = np.load('data/x_coord.npy')
y = np.load('data/y_coord.npy')
z = np.load('data/z_coord.npy')
t = np.load('data/t_coord.npy')

slice = 1
xs_0 = x[0,::slice]                              # x-coordinate of O               
ys_0 = y[0,::slice]                              # y-coordinate of O            
xs_1 = x[1,::slice]                              # x-coordinate of W           
ys_1 = y[1,::slice]                              # y-coordinate of W                         
ts = t[::slice]

dxAB = np.zeros((N_body, nstep), dtype = float)
dyAB = np.zeros((N_body, nstep), dtype = float)
dAB = np.zeros((2,len(xs_0)), dtype = float)

dAB[0,:] = xs_0 - xs_1            
dAB[1,:] = ys_0 - ys_1            
dxAB = x[0,::] - x[1,::]         
dyAB = y[0,::] - y[1,::]          

np.save('data/dxAB', dxAB)     # x-coordinate of binary ow        
np.save('data/dyAB', dyAB)     # y coordinate of binary ow
np.save('data/dAB', dAB)       # (x,y) coordinate of binary ow


