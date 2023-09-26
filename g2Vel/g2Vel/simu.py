import numpy as np
from scipy.interpolate import RectBivariateSpline
from multiprocessing import Pool
from dynesty import DynamicNestedSampler
import pickle
import obstime
import starpara
import rotation
import model
import data


# path to save file
drk = 'data/intfer/465/'

R_a = starpara.R_a
R_b = starpara.R_b
R1 = starpara.R1                  # roche lobe parameter
t = starpara.t
zs = starpara.zs

gX = rotation.gX
gY = rotation.gY
xt = rotation.xt
yt = rotation.yt

sigma = data.sigma1               # sigma1 for lam1 otherwise sigma
sdata = data.sdata1               # sdata1 for lam1 otherwise sdata

lam = obstime.lam
lam1 = obstime.lam1

# orbit fitting coefficients for X and Y
x0, x1, x2, x3 = zs[0], zs[2], zs[4], zs[6]
y0, y1, y2, y3 = zs[1], zs[3], zs[5], zs[7]

# the prior transformation
ini_a, ini_b, ini_1 = R_a, R_b, R1
ini_x0, ini_x1, ini_x2, ini_x3 = x0, x1, x2, x3
ini_y0, ini_y1, ini_y2, ini_y3 = y0, y1, y2, y3
def trans_prior(theta):
    R_ap, R_bp, R1_p, x0_p, x1_p, x2_p, x3_p, y0_p, y1_p, y2_p, y3_p = theta   # include R1_p if observation is for lam1
    
    R_amin, R_amax = 0.5 * ini_a, 1.5 * ini_a 
    R_bmin, R_bmax = 0.5 * ini_b, 1.5 * ini_b
    R1_min, R1_max = 0.5 * ini_1, 1.5 * ini_1                                  # the roche lobe parameter
    
    x0_min, x0_max = 0.5 * ini_x0, 1.5 * ini_x0
    x1_min, x1_max = 0.5 * ini_x1, 1.5 * ini_x1
    x2_min, x2_max = 0.5 * ini_x2, 1.5 * ini_x2
    x3_min, x3_max = 0.5 * ini_x3, 1.5 * ini_x3
    
    y0_min, y0_max = 0.5 * ini_y0, 1.5 * ini_y0
    y1_min, y1_max = 0.5 * ini_y1, 1.5 * ini_y1
    y2_min, y2_max = 0.5 * ini_y2, 1.5 * ini_y2
    y3_min, y3_max = 0.5 * ini_y3, 1.5 * ini_y3
    
    R_a = R_ap * (R_amax - R_amin) + R_amin
    R_b = R_bp * (R_bmax - R_bmin) + R_bmin
    R1 = R1_p * (R1_max - R1_min) + R1_min                                      # the roche lobe parameter
    
    x0 = x0_p * (x0_max - x0_min) + x0_min
    x1 = x1_p * (x1_max - x1_min) + x1_min
    x2 = x2_p * (x2_max - x2_min) + x2_min
    x3 = x3_p * (x3_max - x3_min) + x3_min
    
    y0 = y0_p * (y0_max - y0_min) + y0_min
    y1 = y1_p * (y1_max - y1_min) + y1_min
    y2 = y2_p * (y2_max - y2_min) + y2_min
    y3 = y3_p * (y3_max - y3_min) + y3_min
    
    return (R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3)

# the liklihood function 
def lnlikli(ln):
    global R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3
    R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3 = ln

    g = model.hbt(gX, gY, t, lam1, lam1, R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3)     # use lam=lam1 for lam1 observation                   
    xgr = gX[0,:]                               
    ygr = gY[:,0]                               
    zi = np.transpose(abs(g))                   
    spline = RectBivariateSpline(xgr,ygr,zi)    
    mdata = spline.ev(xt, yt)                   

    G = np.sum(sigma**(-2) * sdata * mdata)                                 
    W = np.sum(sigma**(-2) * mdata * mdata) 

    return 0.5 * (G*G/W - np.log(W))

if __name__ == "__main__":
    
    # the parameters name and number of dimension for parameters
    pnames = ('R_a (ls)', 'R_b (ls)', 'R1 (ls)', 'x0', 'x1', 'x2', 'x3', 'y0', 'y1', 'y2', 'y3')
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
