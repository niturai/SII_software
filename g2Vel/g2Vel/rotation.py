import numpy as np
from obspos import lon, lat, ra, dec, rad, x, y, z
from obstime import jd

# definition of the change in the baseline due to earth rotation 
def rotx(x,y,z,a):
    cs,sn = np.cos(a), np.sin(a)
    return x, cs*y - sn*z, sn*y + cs*z

def roty(x,y,z,a):
    cs,sn = np.cos(a), np.sin(a)
    return cs*x + sn*z, y, -sn*x + cs*z

# Rotate the baseline according to each hour angle
def rotate(dx,dy,dz,jd):
    # Define Hour Angle                                         
    gsid = 18.697374558 + 24.06570982441908*(jd - 2451545)   
    sid = (gsid % 24)*np.pi/12 + lon                          # in 1 hour 15 degree of rotation of earth
    ha = sid - ra
    dx,dy,dz = rotx(dx,dy,dz,-lat)
    dx,dy,dz = roty(dx,dy,dz,ha)
    dx,dy,dz = rotx(dx,dy,dz,dec)
    return dx, dy, dz

# the grids for taking the observation 
def grids(x,y,rad,step):
    N = step
    xup, xdn = np.max(x), np.min(x)                           # xup and xdn will be same as x is a single quantity for a baseline
    xmid, xh = (xup + xdn)/2, (xup - xdn)/2                   # hight (amplitude) of x and mid point
    yup, ydn = np.max(y), np.min(y)
    ymid, yh = (yup + ydn)/2, (yup - ydn)/2                   # hight (amplitude) of y and mid point
    hr = max(xh,yh)                                           # one of the telescope position for a baseline
    r = np.linspace(-hr-2*rad, hr+2*rad, N)                   # (-diameter, diameter) N points on aperture
    wX, wY = np.meshgrid(r,r)                                 # grids of aperture
    gX, gY = xmid + wX, ymid + wY                             # grids for baseline
    return gX, gY


# position of baseline will vary as earth rotate
xt, yt, zt = rotate(x,y,z,jd)                                 # baselines in east and north direction as earth rotate

# grids for each baseline as earth rotate
gX, gY = grids(xt,yt,rad,len(jd))

