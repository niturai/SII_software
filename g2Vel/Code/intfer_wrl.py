import numpy as np
import scipy.special as sp
from scipy.interpolate import RectBivariateSpline
from multiprocessing import Pool
from dynesty import DynamicNestedSampler
import pickle
from obslen import Night, jd, base, ra, dec, lat, lon, lam, lam1, delt, chnl, circ
import math

# the general constant value
c = 299792458.0                  # in m/s
h = 6.62607015e-34               # J.s
k = 1.380649e-23                 # J/K


# the constant parameter of a star
T_a = 35000                      # Temperature of source A
T_b = 57000                        # Temperature of disk
T1 = 20000                      # Temperature of atmosphere
R_a = 39.4503                      # Radius of source A in light seconds
R_b = 13.9236                      # Radius of disk
D = 2.872e10                       # Distance from earth in light seconds

# position of orbit in sky with time
dAB = np.load("data/dAB.npy")
X = dAB[0,0]                  # X-position of orbit
Y = dAB[1,0]                  # Y-position of orbit


# define the radius of emission region
def emis(M1,M2):
    A = np.sqrt(X**2 + Y**2)      # orbital separation of binary
    q = M1/M2                     # mass ratio (M1 having the roche lobe)
    a1 = 0.49*q**(2/3)
    a2 = 0.6*q**(2/3)
    a2 += math.log(1+q**(1/3))
    r1 = A*a1/a2
    return r1


R1 = emis(9,28) - R_b                           # Radius of WR (disk + atmosphere)
print(R1)


# observation Julian day 
jd = jd                           # have (121 shape) 6 minutes apart for Night of observation


# three dimensional baseline 
x = np.real(base)
y = np.imag(base)
z = 1e-6

# radius of aperture
rad = 6     # in meter


# the rotation of x and y baseline due to earth rotation 
def rotx(x,y,z,a):
    cs,sn = np.cos(a), np.sin(a)
    return x, cs*y - sn*z, sn*y + cs*z

def roty(x,y,z,a):
    cs,sn = np.cos(a), np.sin(a)
    return cs*x + sn*z, y, -sn*x + cs*z


# photon flux from single source
def pho_spec(R, T, lam):
    a1 = np.pi * R**2
    a1 /= (D*lam)**2
    a2 = h*c/(lam*k*(T + 2.725))
    a2 = np.clip(a2,0,10)
    I = a1 / (np.exp(a2) - 1)
    return I

# the model 
def hbt(x,y,X,Y,lam,lam1):

    I_a = pho_spec(R_a, T_a, lam)
    I_b = pho_spec(R_b, T_b, lam)
 
    if lam1 == lam:
       I_1 = pho_spec(R1, T1, lam1)
       I_b1 = pho_spec(R_b, T1, lam1)
    else:
       I_1 = I_b1 = 0

    r = np.sqrt(x**2 + y**2)
    v = 2*np.pi/(lam*D)
    
    b1 = v*r*R_a
    b2 = v*r*R_b
    b3 = v*r*R1
    b4 = v*r*R_b

    V_a = 2 * I_a * sp.jv(1, b1)/b1

    V_2 = I_b * sp.jv(1, b2)/b2
    V_3 = I_1 * sp.jv(1, b3)/b3
    V_4 = I_b1 * sp.jv(1, b4)/b4

    V_b = 2 * (V_2 + V_3 - V_4)
    
    V_ab = V_a * V_b * np.cos(v * (x*X + y*Y))

    V = (V_a**2 + V_b**2 + 2 * V_ab)
    V /= (I_a + I_b + I_1 - I_b1)**2
    
    return V


# Rotate the baseline according to each hour angle
def rotate(dx,dy,dz,jd):
    # Define Hour Angle                                         
    gsid = 18.697374558 + 24.06570982441908*(jd - 2451545)   
    sid = (gsid % 24)*np.pi/12 + lon        # in 1 hour 15 degree of rotation of earth
    ha = sid - ra
    dx,dy,dz = rotx(dx,dy,dz,-lat)
    dx,dy,dz = roty(dx,dy,dz,ha)
    dx,dy,dz = rotx(dx,dy,dz,dec)
    return dx,dy


# the grids for taking the observation 
def grids(x,y,rad):
    N = 1024
    xup, xdn = np.max(x), np.min(x)           # xup and xdn will be same as x is a single quantity for a baseline
    xmid, xh = (xup + xdn)/2, (xup - xdn)/2   # hight (amplitude) of x and mid point
    yup, ydn = np.max(y), np.min(y)
    ymid, yh = (yup + ydn)/2, (yup - ydn)/2   # hight (amplitude) of y and mid point
    hr = max(xh,yh)                           # one of the telescope position for a baseline
    r = np.linspace(-hr-2*rad, hr+2*rad, N)   # (-diameter, diameter) N points on aperture
    wX, wY = np.meshgrid(r,r)                 # grids of aperture
    gX, gY = xmid + wX, ymid + wY             # grids for baseline
    return gX, gY, wX, wY


# position of baseline will vary as earth rotate
dist = []
for i in range(len(base)):
    dist.append(rotate(x[i],y[i],z,jd))

distance = np.asarray(dist)

xt = distance[:,0]                            # baselines in east direction as earth rotate
yt = distance[:,1]                            # baselines in north direction as earth rotate


# grids for each baseline as earth rotate
grid = []
for i in range(len(base)):
    grid.append(grids(xt[i,:],yt[i,:],rad))

tgrid = np.asarray(grid)

gX = tgrid[:,0]                                # X axis grids of every baselines

gY = tgrid[:,1]                                # Y axis grids of every baselines

wX = tgrid[:,2]                                # X axis grids of apertures

wY = tgrid[:,3]                                # Y axis grids of apertures

# The aperture function 
w = circ(wX,wY,rad)

'''
# the total flux spectrum from source
def tot_spec(lam,lam1):
    sp = pho_spec(R_a, T_a, lam)
    sp += pho_spec(R_a, T_b, lam)
    if lam1 == lam:
       sp += pho_spec(R1, T1, lam1)
       sp -= pho_spec(R_b, T1, lam1)
    else:
       sp += 0
    return sp
    
tspec = tot_spec(lam, lam1)       # for continuous spectra
tspec1 = tot_spec(lam1, lam1)     # for line spectra


print('total flux spectrum from continuous line=', tspec)
print('total flux spectrum from continuous line=', tspec1)
'''
tspec = 4.9e-5                         # for continuous spectra
tspec1 = 0.98e-4                       # for line spectra


# The value of area of aperture
area = np.pi * rad**2
print("area =",area)


# the noise from each telescope with 1 sec of average time of observation
sigma = (area * tspec)**(-1) * np.sqrt(delt)/chnl                         # for continuous spectra
sigma1 = (area * tspec1)**(-1) * np.sqrt(delt)/chnl                       # for line spectra
print('noise =',sigma)
print('noise =',sigma1)


# simulate the signal for all baseline
def smdata(sigma,lam,lam1):
    
    g = hbt(gX,gY,X,Y,lam,lam1)               # signal for all baselines,  for disk only
    
    # simulate the data
    simu_data = []
    
    for i in range(len(base)):
        gi = g[i,:,:]                             # signal for ith baseline
        xgr = gX[i,0,:]                           # 0th elements of each x axis
        ygr = gY[i,:,0]                           # 0th elements of each y axis
        zi = np.transpose(abs(gi))                # take a transpose of absolute value of signal so that it will belongs to each (xgr, ygr)
        spline = RectBivariateSpline(xgr,ygr,zi)  # interpolation and smoothing of data zi
        simu_data.append(spline.ev(xt[i,:], yt[i,:]) + sigma * np.random.randn(len(jd)))                        # Returns the interpolated value for baseline position (xt[i,:], yt[i,:])
    
    return np.array(simu_data)     

sdata = smdata(sigma,lam,lam1)
sdata1 = smdata(sigma1,lam1,lam1)


# total SNR value
SNR = (1/sigma)*(np.sum(sdata)**2)**.5
print("total signal to noise ratio =", SNR)


if __name__ == "__main__":
     # save the signal and observation time
     np.save("data/jd",jd)
     np.save("data/xt",xt)
     np.save("data/yt",yt)
     np.save("data/sig",smdata(0,lam,lam1))
     np.save("data/sig_noise",smdata(sigma,lam,lam1))
     np.save("data/sig1",smdata(0,lam1,lam1))
     np.save("data/sig_noise1",smdata(sigma,lam1,lam1))

# the prior transformation
ini_x, ini_y = R_a, R_b
ini_X, ini_Y = X, Y
def trans_prior(theta):
    R_aprime, R_bprime, X_prime, Y_prime = theta
    
    R_amin, R_amax = 0.5 * ini_x, 1.5 * ini_x
      
    R_bmin, R_bmax = 0.5 * ini_y, 1.5 * ini_y


    X_min, X_max = 0.5 * ini_X, 1.5 * ini_X

    Y_min, Y_max = 0.5 * ini_Y, 1.5 * ini_Y
     
    
    R_a = R_aprime * (R_amax - R_amin) + R_amin
    R_b = R_bprime * (R_bmax - R_bmin) + R_bmin

    X = X_prime * (X_max - X_min) + X_min
    Y = Y_prime * (Y_max - Y_min) + Y_min
    
    return (R_a, R_b, X, Y)

# the liklihood function 
def lnlikli(ln):
    global R_a, R_b, X, Y
    R_a, R_b, X, Y = ln
    g = hbt(gX,gY,X,Y,lam,lam1)                # Signal for all baseline
    model_data = []
    for i in range(len(base)):
        gi = g[i,:,:]
        xgr = gX[i,0,:]
        ygr = gY[i,:,0]
        zi = np.transpose(abs(gi))
        spline = RectBivariateSpline(xgr,ygr,zi)
        model_data.append(spline.ev(xt[i,:],yt[i,:])) 
        
    mdata = np.array(model_data)

    G = np.sum(sigma**(-2) * sdata * mdata)
    W = np.sum(sigma**(-2) * mdata * mdata) 

    return 0.5 * (G*G/W - np.log(W))

if __name__ == "__main__": 
    # the parameters name and number of dimension for parameters
    pnames = ('$R_a (ls) $', '$R_b (ls) $', 'X', 'Y')
    ndim = len(pnames)

    # the sampling function 
    pool = Pool()
    sampler = DynamicNestedSampler(lnlikli, trans_prior, ndim, pool=pool, queue_size=28)
    sampler.run_nested()

    # Save the sample data in form of pickle 
    dres = sampler.results
    out = open('data/dres_wrl.pkl','wb')
    pickle.dump(dres,out)
    out.close()
    pool.close()









