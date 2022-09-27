import numpy as np
import scipy.special as sp
from scipy.interpolate import RectBivariateSpline
from multiprocessing import Pool
from dynesty import DynamicNestedSampler
import pickle
from obslen import Night, jd, base, ra, dec, lat, lon, lam, delt, chnl, circ


# the general constant value
c = 299792458.0                  # in m/s
h = 6.62607015e-34               # J.s
k = 1.380649e-23                 # J/K


# the constant parameter of Spica binary system
T_a = 35000                        # Temperature of source A in Kelvin
T_b = 57000                        # Temperature of source B in Kelvin
R_a = 39.4503                      # Radius of source A in light seconds
R_b = 13.9236                      # Radius of source B in light seconds
D = 2.872e10                       # Distance from earth in light seconds 


# position of orbit in sky with time
dAB = np.load("data/dAB.npy")
X = dAB[0,0]                  # X-position of orbit
Y = dAB[1,0]                  # Y-position of orbit

# limb-darkening coefficients
l = 0.6
m = 0.2
n = 0.1

# parameters, which are suppose to be estimated
truth = [R_a,R_b, X, Y]
print("Parameters of both stars =",truth)

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


# photon flux from source A or B
def pho_spec(R, T, lam):
    a1 = np.pi * (R/D)**2
    a1 /= lam**2
    a2 = h*c/(lam*k*(T + 2.725))
    a2 = np.clip(a2,0,10)
    I = a1 / (np.exp(a2) - 1) 
    return I

# the model (which is HBT correlation) to estimate the parameter
def hbt(x,y,X,Y,lam):
    r = np.sqrt(x**2 + y**2)
    v = 2*np.pi/(lam*D)

    I_a = pho_spec(R_a,T_a,lam)
    n_a = v*r*R_a

    I_b = pho_spec(R_b,T_b,lam)
    n_b = v*r*R_b   

    # function of limb-darkening coefficient
    C = 1/2 - l/6 - m/12 - n/20
    a = 1- l - m - n
    b = l + 2*m + 3*n
    f = -m - 3*n
    g = n

    t1 = a*sp.jv(1,n_a)/n_a
    t2 = b*np.sqrt(np.pi/2)*sp.jv(3/2,n_a)/n_a**(3/2)
    t3 = 2*f*sp.jv(2,n_a)/n_a**2
    t4 = 3*g*np.sqrt(np.pi/2)*sp.jv(5/2,n_a)/n_a**(5/2)
    t = t1 + t2 + t3 + t4
 
    s1 = a*sp.jv(1,n_b)/n_b
    s2 = np.sqrt(np.pi/2)*b*sp.jv(3/2,n_b)/n_b**(3/2)
    s3 = 2*f*sp.jv(2,n_b)/n_b**2
    s4 = 3*g*np.sqrt(np.pi/2)*sp.jv(5/2,n_b)/n_b**(5/2)
    s = s1 + s2 + s3 + s4

    V_a = I_a * t  
    V_b = I_b * s
    V_ab = V_a * V_b * np.cos(v * (x*X + y*Y))

    return (I_a*C + I_b*C)**(-2) * (V_a**2 + V_b**2 + 2*V_ab) 

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


# the total flux spectrum from source
tspec = pho_spec(R_a,T_a,lam) + pho_spec(R_a, T_b,lam)
print('total flux spectrum =', tspec)


# The value of area of aperture
area = np.pi * rad**2
print("area =",area)


# the noise from each telescope
sigma = (area * tspec)**(-1) * np.sqrt(delt)/chnl
print('noise =',sigma)


# the signal without any noise
g = hbt(gX,gY,X,Y,lam)                            # signal for all baselines


# simulate the signal for all baseline
def smdata(sigma):
    simu_data = []
    for i in range(len(base)):
        gi = g[i,:,:]                             # signal for ith baseline
        xgr = gX[i,0,:]                           # 0th elements of each x axis
        ygr = gY[i,:,0]                           # 0th elements of each y axis
        zi = np.transpose(abs(gi))                # take a transpose of absolute value of signal so that it will belongs to each (xgr, ygr)
        spline = RectBivariateSpline(xgr,ygr,zi)  # interpolation and smoothing of data zi
        simu_data.append(spline.ev(xt[i,:], yt[i,:]) + sigma * np.random.randn(len(jd)))                        # Returns the interpolated value for baseline position (xt[i,:], yt[i,:])
    
    return np.array(simu_data)     

sdata = smdata(sigma)

# save the signal and observation time
np.save("data/jd",jd)
np.save("data/xt",xt)
np.save("data/yt",yt)
np.save("data/sig",smdata(0))
np.save("data/sig_noise",smdata(sigma))


# total SNR value
SNR = (1/sigma)*(np.sum(sdata)**2)**.5
print("total signal to noise ratio =", SNR)


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
    g = hbt(gX,gY,X,Y,lam)                # Signal for all baseline
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
    out = open('data/dres.pkl','wb')
    pickle.dump(dres,out)
    out.close()
    pool.close()

