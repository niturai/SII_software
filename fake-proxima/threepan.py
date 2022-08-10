from numpy import hstack
import imageio

def run(lp,rp,F):
    images = []
    lim = imageio.imread('tmp/src.png')
    mim = imageio.imread('tmp/corr.png')
    for f in range(F):
        rname = 'tmp/cdiff%i.png' % f
        rim = imageio.imread(rname)
        print('read image ',f);
        im = hstack((lim,mim,rim))
        images.append(im)
    imageio.mimsave('chzcen-hess.gif', images, fps=5)

    
run('sky','cdens',76)


