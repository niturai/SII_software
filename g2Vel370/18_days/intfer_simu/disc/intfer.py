import numpy as np
import scipy.special as sp
from scipy.interpolate import RectBivariateSpline
from multiprocessing import Pool
from dynesty import DynamicNestedSampler
import pickle
from obslen import Night, jd, base, ra, dec, lat, lon, lam, lam1, delt
import math

# the general constant value
c = 299792458.0                  # in m/s
h = 6.62607015e-34               # J.s
k = 1.380649e-23                 # J/K


# the constant parameter of a star
T_a = 35000                      # Temperature of source A
T_b = 57000                        # Temperature of disk
R_a = 39.4503                      # Radius of source A in light seconds
R_b = 13.9236                      # Radius of disk
D = 3.458e10                       # Distance from earth in light seconds

# position of orbit in sky with time
dAB = np.load("../../orbt_simu/data/dAB.npy")
Xth = []
Yth = []
for i in range(0, Night, 1):
    Xth.append(dAB[0, (21+27)*i : 21*(i+1)+27*i])
    Yth.append(dAB[1, (21+27)*i : 21*(i+1)+27*i])
    

X = np.array(Xth)                           # X-position of orbit for 18 days and 21 steps of observation b/w 7pm to 5am
Y = np.array(Yth)                           # Y-position of orbit for 18 days and 21 steps of observation b/w 7pm to 5am

# three dimensional baseline 
x = round(np.real(base),2)
y = round(np.imag(base),2)
z = 1e-6


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
def hbt(x,y,X,Y,lam):

    I_a = pho_spec(R_a, T_a, lam)
    I_b = pho_spec(R_b, T_b, lam)

    r = np.sqrt(x**2 + y**2)
    v = 2*np.pi/(lam*D)
    
    b1 = v*r*R_a
    b2 = v*r*R_b

    V_a = 2 * I_a * sp.jv(1, b1)/b1

    V_b = 2 * I_b * sp.jv(1, b2)/b2
    
    V_ab = V_a * V_b * np.cos(v * (x*X + y*Y))

    V = (V_a**2 + V_b**2 + 2 * V_ab)
    V /= (I_a + I_b)**2
    
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

rad = 6
# the grids for taking the observation 
def grids(x,y,rad):
    N = 21
    xup, xdn = np.max(x), np.min(x)           # xup and xdn will be same as x is a single quantity for a baseline
    xmid, xh = (xup + xdn)/2, (xup - xdn)/2   # hight (amplitude) of x and mid point
    yup, ydn = np.max(y), np.min(y)
    ymid, yh = (yup + ydn)/2, (yup - ydn)/2   # hight (amplitude) of y and mid point
    hr = max(xh,yh)                           # one of the telescope position for a baseline
    r = np.linspace(-hr-2*rad, hr+2*rad, N)   # (-diameter, diameter) N points on aperture
    wX, wY = np.meshgrid(r,r)                 # grids of aperture
    gX, gY = xmid + wX, ymid + wY             # grids for baseline
    return gX, gY


# position of baseline will vary as earth rotate
dist = []
for i in range(0, Night, 1):
    dist.append(rotate(x,y,z,jd[i]))

distance =np.asarray(dist)          

xt = distance[:,0]                  # baselines in east direction as earth rotate
yt = distance[:,1]                  # baselines in north direction as earth rotate


# grids for each baseline as earth rotate
grid = []                           
for i in range(0, Night, 1):
    grid.append(grids(xt[i], yt[i], rad))
    
tgrid = np.asarray(grid)

gX = tgrid[:,0]                    # grids in X axix for 18 days

gY = tgrid[:,1]                    # grids in Y axix for 18 days


# total photon spectrun in effective area 
tspec = 2e-5                                             # for continuous spectra A*\Phi

# the noise from each telescope with 30 min of average time of observation
sigma = np.sqrt(delt/1800)/tspec                         

print('noise =',sigma)


# simulate the signal for all baseline
def smdata(sigma,lam):
    sm_data = []
    for i in range(0, Night, 1):
        gi = hbt(gX[i,:],gY[i,:],X[i],Y[i],lam)                       
        xgr = gX[i,0,:]                           # 0th elements of each x axis
        ygr = gY[i,:,0]                           # 0th elements of each y axis
        zi = np.transpose(abs(gi))                # take a transpose of absolute value of signal so that it will belongs to each (xgr, ygr)
        spline = RectBivariateSpline(xgr,ygr,zi)  # interpolation and smoothing of data zi
        sm_data.append(spline.ev(xt[i], yt[i]) + sigma * np.random.randn(len(jd[i])))    # Returns the interpolated value for baseline position (xt,yt)
    
    return np.array(sm_data)     

sdata = smdata(sigma,lam)

# total SNR value
SNR = (1/sigma)*np.sum(sdata**2)**.5
print("total signal to noise ratio =", SNR)


# the prior transformation
ini_x, ini_y = R_a, R_b
def trans_prior(theta):
    R_aprime, R_bprime = theta
    
    R_amin, R_amax = 0.5 * ini_x, 1.5 * ini_x
      
    R_bmin, R_bmax = 0.5 * ini_y, 1.5 * ini_y

    
    R_a = R_aprime * (R_amax - R_amin) + R_amin
    R_b = R_bprime * (R_bmax - R_bmin) + R_bmin
    
    return (R_a, R_b)

# the liklihood function 
def lnlikli(ln):
    global R_a, R_b
    R_a, R_b = ln
    m_data = []
    for i in range(0, Night, 1):
        gi = hbt(gX[i,:],gY[i,:],X[i],Y[i],lam)                               # for ith day
        xgr = gX[i,0,:]                           # 0th elements of each x axis
        ygr = gY[i,:,0]                           # 0th elements of each y axis
        zi = np.transpose(abs(gi))                # take a transpose of absolute value of signal so that it will belongs to each (xgr, ygr)
        spline = RectBivariateSpline(xgr,ygr,zi)  # interpolation and smoothing of data zi
        m_data.append(spline.ev(xt[i], yt[i]))          # Returns the interpolated value for baseline position (xt,yt)
    mdata = np.array(m_data)
    G = np.sum(sigma**(-2) * sdata * mdata)
    W = np.sum(sigma**(-2) * mdata * mdata) 

    return 0.5 * (G*G/W - np.log(W))

if __name__ == "__main__": 
    # save the signal and observation time
    np.save("data/jd",jd)
    np.save("data/xt",xt)
    np.save("data/yt",yt)
    np.save("data/sig",smdata(0,lam))
    np.save("data/sig_noise",smdata(sigma,lam))

    # the parameters name and number of dimension for parameters
    pnames = ('$R_a (ls) $', '$R_b (ls) $')
    ndim = len(pnames)

    # the sampling function 
    pool = Pool()
    sampler = DynamicNestedSampler(lnlikli, trans_prior, ndim, pool=pool, queue_size=28)
    sampler.run_nested()

    # Save the sample data in form of pickle 
    dres = sampler.results
    out = open('data/dres.pkl','wb')
    pickle.dump(dres,out)
    out.close()
    pool.close()

