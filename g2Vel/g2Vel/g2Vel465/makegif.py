import imageio.v2 as imageio
from julian import Night

# number of png files in each folder
N = Night

# path to read the png files
dfig = ['fig/sig/sig', 'fig/cntr/cntr']

# path to save gif movie files 
dgif = ['gif/sig', 'gif/cntr']


for i in range(0, len(dfig), 1):
    images = []
    for f in range(1,N+1,1):
        fname = dfig[i] + ('%i' % f) + '.png'    
        images.append(imageio.imread(fname))
    imageio.mimsave(dgif[i] + '.gif', images, fps=1)

