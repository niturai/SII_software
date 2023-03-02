import imageio.v2 as imageio

def run(ops,F):
    images = []
    for f in range(1,F,1):
        fname = 'fig/370/cntr/'+ops+('%i' % f)+'.png'    
        images.append(imageio.imread(fname))
    imageio.mimsave('gif/' + ops + '370' +'.gif', images, fps=1)


run('cntr',17)


