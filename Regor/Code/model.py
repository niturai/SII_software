import numpy as np
import math
import scipy.special as sp
from scipy.interpolate import RectBivariateSpline
from multiprocessing import Pool
from dynesty import DynamicNestedSampler
import pickle
from obspos import A
from obstime import lam, lam1, ita, loss, delt, jd
from julian import dt
from starpara import c, h, k, T_a, T_b, T1, R_a, R_b, R1, D, X, Y
from rotation import xt, yt, gX, gY

# photon flux from single source
def pho_spec(R, T, lam):
    a1 = np.pi * R**2
    a1 /= (D*lam)**2
    a2 = h*c/(lam*k*(T + 2.725))
    a2 = np.clip(a2,0,10)
    I = a1 / (np.exp(a2) - 1)
    return I


# define the model of signal
def hbt(x,y,X,Y,lam,lam1):                       # lam is for continuum and lam1 is for line spectrum

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

# define the total photon spectrum for lam (A*\Phi)
def tot_spec(lam,lam1):
    sp = pho_spec(R_a, T_a, lam)
    sp += pho_spec(R_b, T_b, lam)
    if lam1 == lam:
       sp += pho_spec(R1, T1, lam1)
       sp -= pho_spec(R_b, T1, lam1)
    else:
       sp += 0
    return sp
    
tspec = A*ita*loss*tot_spec(lam, lam1)       # for lam spectra
tspec1 = 2e-5                                  # for lam1 spectra    

# total flux 
print('total flux for 370nm',tspec)
print('total flux for 465nm',tspec1)


# the noise for each baseline
sigma = np.sqrt(delt/dt)/tspec            # for lam  
sigma1 = np.sqrt(delt/dt)/tspec1          # for lam1
print('noise for 370nm line =',sigma)
print('noise for 465nm line =',sigma1)

# simulate the signal for all baseline
def smdata(sigma,lam,lam1):
    # the signal for all days
    g = hbt(gX,gY,X,Y,lam,lam1)                       
    xgr = gX[0,:]                           # 0th elements of each x axis
    ygr = gY[:,0]                           # 0th elements of each y axis
    zi = np.transpose(abs(g))                # take a transpose of absolute value of signal so that it will belongs to each (xgr, ygr)
    spline = RectBivariateSpline(xgr,ygr,zi)  # interpolation and smoothing of data zi
    data = spline.ev(xt, yt) + sigma * np.random.randn(len(jd)) # Returns the interpolated value for baseline position (xt,yt)
    
    return data     

sdata = smdata(sigma,lam,lam1)
sdata1 = smdata(sigma1,lam1,lam1)

# calculate the Signal-to-noise ratio
SNR = sigma**(-1)*np.sum(sdata**2)**.5
SNR1 = sigma1**(-1)*np.sum(sdata1**2)**.5
print("total signal to noise ratio for 370nm =", SNR)
print("total signal to noise ratio for 465nm =", SNR1)

# the prior transformation
ini_x, ini_y, ini_1 = R_a, R_b, R1
def trans_prior(theta):
    R_aprime, R_bprime, R_1prime = theta
    
    R_amin, R_amax = 0.5 * ini_x, 1.5 * ini_x     
    R_bmin, R_bmax = 0.5 * ini_y, 1.5 * ini_y
    R_1min, R_1max = 0.5 * ini_1, 1.5 * ini_1
    

    
    R_a = R_aprime * (R_amax - R_amin) + R_amin
    R_b = R_bprime * (R_bmax - R_bmin) + R_bmin
    R1 = R_1prime * (R_1max - R_1min) + R_1min
    
    return (R_a, R_b, R1)

# the liklihood function 
def lnlikli(ln):
    global R_a, R_b, R1
    R_a, R_b, R1 = ln

    # the signal for all days
    g = hbt(gX,gY,X,Y,lam1,lam1)                       
    xgr = gX[0,:]                           # 0th elements of each x axis
    ygr = gY[:,0]                           # 0th elements of each y axis
    zi = np.transpose(abs(g))                # take a transpose of absolute value of signal so that it will belongs to each (xgr, ygr)
    spline = RectBivariateSpline(xgr,ygr,zi)  # interpolation and smoothing of data zi
    mdata = spline.ev(xt, yt)               # Returns the interpolated value for baseline position (xt,yt)

    G = np.sum(sigma1**(-2) * sdata1 * mdata)
    W = np.sum(sigma1**(-2) * mdata * mdata) 

    return 0.5 * (G*G/W - np.log(W))


# save the data file
if __name__ == "__main__":
   # save the signal and observation time
   np.save("data/465/jd",jd)
   np.save("data/465/xt",xt)
   np.save("data/465/yt",yt)
   np.save("data/465/sig",smdata(0,lam1,lam1))
   np.save("data/465/sig_noise",smdata(sigma1,lam1,lam1))


  # the parameters name and number of dimension for parameters
   pnames = ('$R_a (ls)$', '$R_b (ls)$', '$R1 (ls)$')
   ndim = len(pnames)

 # the sampling function 
   pool = Pool()
   sampler = DynamicNestedSampler(lnlikli, trans_prior, ndim, pool=pool, queue_size=28)
   sampler.run_nested()

 # Save the sample data in form of pickle 
   dres = sampler.results
   out = open('data/465/dres2.pkl','wb')
   pickle.dump(dres,out)
   out.close()
   pool.close()

