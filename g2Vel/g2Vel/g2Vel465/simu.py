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
drk = 'data/'

R_a = starpara.R_a
R_b = starpara.R_b
R1 = starpara.R1                  # roche lobe parameter
X = starpara.X
Y = starpara.Y

gX = rotation.gX
gY = rotation.gY
xt = rotation.xt
yt = rotation.yt

sigma = data.sigma1               # sigma1 for lam1 otherwise sigma
sdata = data.sdata1               # sdata1 for lam1 otherwise sdata

lam = obstime.lam
lam1 = obstime.lam1

# the prior transformation
ini_a, ini_b, ini_1 = R_a, R_b, R1
def trans_prior(theta):
    R_ap, R_bp, R1_p = theta   # include R1_p if observation is for lam1
    
    R_amin, R_amax = 0.5 * ini_a, 1.5 * ini_a 
    R_bmin, R_bmax = 0.5 * ini_b, 1.5 * ini_b
    R1_min, R1_max = 0.5 * ini_1, 1.5 * ini_1                                  # the roche lobe parameter
    
    R_a = R_ap * (R_amax - R_amin) + R_amin
    R_b = R_bp * (R_bmax - R_bmin) + R_bmin
    R1 = R1_p * (R1_max - R1_min) + R1_min                                     # the roche lobe parameter
    
    return (R_a, R_b, R1)

# the liklihood function 
def lnlikli(ln):
    global R_a, R_b, R1
    R_a, R_b, R1 = ln

    g = model.hbt(gX, gY, X, Y, lam1, lam1, R_a, R_b, R1)                         # use lam=lam1 for lam1 observation                   
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
    pnames = ('R_a (ls)', 'R_b (ls)', 'R1 (ls)')
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
