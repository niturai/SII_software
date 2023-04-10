import numpy as np
import math
import scipy.special as sp
from scipy.interpolate import RectBivariateSpline
from multiprocessing import Pool
from dynesty import DynamicNestedSampler
import pickle
from julian import dt
from obstime import jd, delt, lam, lam1, ita, loss
from obspos import A
from starpara import c, h, k, T_a, T_b, T1, R_a, R_b, R1, D, t, zs
from rotation import xt, yt, gX, gY

# path to save file
drk = 'data/intfer/465/'

# orbit fitting coefficients for X and Y
x0, x1, x2 = zs[0], zs[2], zs[4]
y0, y1, y2 = zs[1], zs[3], zs[5]

# photon flux from single source
def pho_spec(R, T, lam):
    a1 = np.pi * R**2
    a1 /= (D*lam)**2
    a2 = h*c/(lam*k*(T + 2.725))
    a2 = np.clip(a2,0,10)
    I = a1 / (np.exp(a2) - 1)
    return I

# define the model of signal
def hbt(x,y,t,lam,lam1):                       # lam is for continuum and lam1 is for line spectrum

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

    X = x0 + x1*t + x2*t**2
    Y = y0 + y1*t + y2*t**2

    V_a = 2 * I_a * sp.jv(1, b1)/b1

    V_2 = I_b * sp.jv(1, b2)/b2
    V_3 = I_1 * sp.jv(1, b3)/b3
    V_4 = I_b1 * sp.jv(1, b4)/b4

    V_b = 2 * (V_2 + V_3 - V_4)
    
    V_ab = V_a * V_b * np.cos(v * (x*X + y*Y))

    V = (V_a**2 + V_b**2 + 2 * V_ab)
    V /= (I_a + I_b + I_1 - I_b1)**2
    
    return V

# the photon spectrum (\Phi) for source as disk
def tot_spec(lam,lam1):
    sp = pho_spec(R_a, T_a, lam)
    sp += pho_spec(R_b, T_b, lam)
    if lam1 == lam:
       sp += pho_spec(R1, T1, lam1)
       sp -= pho_spec(R_b, T1, lam1)
    else:
       sp += 0
    return sp

# total photon spectrum on HESS observatory
tspec = A*ita*loss*tot_spec(lam, lam1)         # for lam spectra, with 50% loss in sysytem
tspec1 = 2e-5                                  # for lam1 spectra      

# the noise for each baseline
sigma = np.sqrt(delt/dt)/tspec                 # for lam  
sigma1 = np.sqrt(delt/dt)/tspec1               # for lam1
print('noise for 370nm line =',sigma)
print('noise for 465nm line =',sigma1)

# simulate the observed signal for all baseline
def smdata(sigma,lam,lam1):
    g = hbt(gX,gY,t,lam,lam1)                 # the signal for all days

    # convert grids into x and y                       
    xgr = gX[0,:]                             # 0th elements of each y
    ygr = gY[:,0]                             # 0th elements of each x  
    zi = np.transpose(abs(g))                 # Now g belongs to each (xgr, ygr)
    spline = RectBivariateSpline(xgr,ygr,zi)  # interpolation to do smoothing of data zi
    data = spline.ev(xt, yt) + sigma * np.random.randn(len(jd)) # Returns the interpolated value for baseline position (xt,yt)
    return data      

sdata = smdata(sigma,lam,lam1)                    # for continuous line
sdata1 = smdata(sigma1,lam1,lam1)                 # for carbon line

# calculate the Signal-to-noise ratio
SNR = sigma**(-1)*np.sum(sdata**2)**.5
SNR1 = sigma1**(-1)*np.sum(sdata1**2)**.5
print("total signal to noise ratio for 370nm =", SNR)
print("total signal to noise ratio for 465nm =", SNR1)
   

# the prior transformation
ini_a, ini_b, ini_1 = R_a, R_b, R1
ini_x0, ini_x1, ini_x2 = x0, x1, x2
ini_y0, ini_y1, ini_y2 = y0, y1, y2
def trans_prior(theta):
    R_ap, R_bp, R_1p, x0_p, x1_p, x2_p, y0_p, y1_p, y2_p = theta
    
    R_amin, R_amax = 0.5 * ini_a, 1.5 * ini_a 
    R_bmin, R_bmax = 0.5 * ini_b, 1.5 * ini_b
    R_1min, R_1max = 0.5 * ini_1, 1.5 * ini_1
    
    x0_min, x0_max = 0.5 * ini_x0, 1.5 * ini_x0
    x1_min, x1_max = 0.5 * ini_x1, 1.5 * ini_x1
    x2_min, x2_max = 0.5 * ini_x2, 1.5 * ini_x2
    
    y0_min, y0_max = 0.5 * ini_y0, 1.5 * ini_y0
    y1_min, y1_max = 0.5 * ini_y1, 1.5 * ini_y1
    y2_min, y2_max = 0.5 * ini_y2, 1.5 * ini_y2
    
    R_a = R_ap * (R_amax - R_amin) + R_amin
    R_b = R_bp * (R_bmax - R_bmin) + R_bmin
    R1 = R_1p * (R_1max - R_1min) + R_1min
    
    x0 = x0_p * (x0_max - x0_min) + x0_min
    x1 = x1_p * (x1_max - x1_min) + x1_min
    x2 = x2_p * (x2_max - x2_min) + x2_min
    
    y0 = y0_p * (y0_max - y0_min) + y0_min
    y1 = y1_p * (y1_max - y1_min) + y1_min
    y2 = y2_p * (y2_max - y2_min) + y2_min
    
    return (R_a, R_b, R1, x0, x1, x2, y0, y1, y2)

# the liklihood function 
def lnlikli(ln):
    global R_a, R_b, R1, x0, x1, x2, y0, y1, y2
    R_a, R_b, R1, x0, x1, x2, y0, y1, y2 = ln

    # the signal for all days
    g = hbt(gX,gY,t,lam1,lam1)                        
    xgr = gX[0,:]                               
    ygr = gY[:,0]                               
    zi = np.transpose(abs(g))                   
    spline = RectBivariateSpline(xgr,ygr,zi)    
    mdata = spline.ev(xt, yt)                   

    G = np.sum(sigma1**(-2) * sdata1 * mdata)
    W = np.sum(sigma1**(-2) * mdata * mdata) 

    return 0.5 * (G*G/W - np.log(W))


if __name__ == "__main__":
 
    # save the signal and observation time
    np.save(drk + "jd", jd)
    np.save(drk + "xt", xt)
    np.save(drk + "yt", yt)
    np.save(drk + "sig", smdata(0,lam1,lam1))
    np.save(drk + "sig_noise", smdata(sigma1,lam1,lam1))
   
    # the parameters name and number of dimension for parameters
    pnames = ('R_a (ls)', 'R_b (ls)', 'R_b (ls)', 'x0', 'x1', 'x2', 'y0', 'y1', 'y2')
    ndim = len(pnames)

    # the sampling function 
    pool = Pool()
    sampler = DynamicNestedSampler(lnlikli, trans_prior, ndim, pool=pool, queue_size=28)
    sampler.run_nested()

    # Save the sample data in form of pickle 
    dres = sampler.results
    out = open(drk + 'dres.pkl','wb')
    pickle.dump(dres,out)
    out.close()
    pool.close()

