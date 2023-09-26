import numpy as np
from scipy.interpolate import RectBivariateSpline
from julian import dt
from obstime import delt, lam, lam1, ita, loss, jd
from obspos import A
from rotation import xt, yt, zt, gX, gY
from starpara import c, h, k, T_a, T_b, T1, R_a, R_b, R1, D, t, zs
from model import hbt, pho_spec


# path to save file
drk = 'data/intfer/465/'

# orbit fitting coefficients for X and Y
x0, x1, x2, x3 = zs[0], zs[2], zs[4], zs[6]
y0, y1, y2, y3 = zs[1], zs[3], zs[5], zs[7]

# define the total photon spectrum \Phi for lam or lam1
def tot_spec(lam,lam1):
    sp = pho_spec(R_a, T_a, lam)
    sp += pho_spec(R_b, T_b, lam)
    if lam1 == lam:
       sp += pho_spec(R1, T1, lam1)
       sp -= pho_spec(R_b, T1, lam1)
    else:
       sp += 0
    return sp
    
tspec = A*ita*loss*tot_spec(lam, lam1)                         # for lam spectra, with 50% loss in sysytem
tspec1 = 2e-5                                                  # for lam1 spectra      

# the noise for each baseline
sigma = np.sqrt(delt/dt)/tspec                                 # for lam  
sigma1 = np.sqrt(delt/dt)/tspec1                               # for lam1
print('noise for 370nm line =',sigma)
print('noise for 465nm line =',sigma1)

# simulate the signal for all baseline
def smdata(sigma,lam,lam1):
    g = hbt(gX, gY, t, lam, lam1, R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3)                    # the signal for all days      
    xgr = gX[0,:]                                               # 0th elements of each x axis
    ygr = gY[:,0]                                               # 0th elements of each y axis
    zi = np.transpose(abs(g))                                   # it will belongs to each (xgr, ygr)
    spline = RectBivariateSpline(xgr,ygr,zi)                    # interpolation and smoothing of data zi
    data = spline.ev(xt, yt) + sigma * np.random.randn(len(jd)) # the interpolated value for baseline (xt,yt) with gaussian distribution noise 
    return data     

sdata = smdata(sigma,lam,lam1)                                  # for continuous line
sdata1 = smdata(sigma1,lam1,lam1)                               # for carbon line

# calculate the Signal-to-noise ratio
SNR = sigma**(-1)*np.sum(sdata**2)**.5
SNR1 = sigma1**(-1)*np.sum(sdata1**2)**.5
print("total signal to noise ratio for 370nm =", SNR)
print("total signal to noise ratio for 465nm =", SNR1)

# save the data file
if __name__ == "__main__":
   # save the signal and observation time
   np.save(drk + "jd", jd)
   np.save(drk + "xt", xt)
   np.save(drk + "yt", yt)
   np.save(drk + "sig", smdata(0,lam1,lam1))                     # model data for observation wavelength
   np.save(drk + "sig_noise", smdata(sigma1,lam1,lam1))           # noisy data for observation wavelength


