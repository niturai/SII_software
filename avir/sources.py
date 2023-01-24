import numpy as np

# Produces an elliptical source
def blob(sx,sy,rad,posx=0,posy=0,elong=1,pa=0):
    dx,dy = sx-posx, sy-posy
    cs,sn = np.cos(pa),-np.sin(pa)
    dx,dy = cs*dx - sn*dy, sn*dx + cs*dy
    bla = 0*sx
    bla[dx**2*elong + dy**2/elong < rad**2] = 1
    return bla

# Produces a crescent source
def crescent(sx,sy,rad,f):
    cres = 0*sx
    cres[sx**2 + sy**2 < rad**2] = 1
    cres[(sx-f*rad)**2 + sy**2 < ((1-f)*rad)**2] = 0
    return cres

# Smooth to reduce pixellation artefacts
def smooth(S,W):
    M = 0*S
    ksum = 0
    for i in range(1-W,W):
        for j in range(1-W,W):
            ker = 1/(W*W + i*i + j*j)
            ksum += ker
            M[W:-W,W:-W] += ker*S[W+i:-W+i,W+j:-W+j]
    M /= ksum
    return M

