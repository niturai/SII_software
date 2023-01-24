from numpy import hstack
import imageio

def run(lp,rp,F):
    images = []
    for f in range(1,F):
        lname = 'tmp/'+lp+('%i' % f)+'.png'
        rname = 'tmp/'+rp+('%i' % f)+'.png'
        lim = imageio.imread(lname)
        rim = imageio.imread(rname)
        print('read image ',f);
        im = hstack((lim,rim))
        #imageio.imsave(fname,im)
        images.append(im)
    imageio.mimsave('spica.gif', images, fps=5)

    
run('sky','cdens',10)


