import numpy as np
from dynesty.utils import resample_equal
import pickle
import matplotlib.pyplot as plt
import corner
import math
from starpara import R_a, R_b, R1, zs

# path to read file
drk = 'data/intfer/465/'

# path to save file
dfig = 'fig/intfer/465/'

# orbit fitting coefficients for X and Y
x0, x1, x2, x3 = zs[0], zs[2], zs[4], zs[6]
y0, y1, y2, y3 = zs[1], zs[3], zs[5], zs[7]

# parameters, which are estimated
truth = [R_a, R_b, R1, x0, x1, x2, x3, y0, y1, y2, y3]
print('parameters',truth)


#number of parameter to estimate
pnames = ('R_a (ls)', 'R_b (ls)', 'R1 (ls)', 'x0', 'x1', 'x2', 'x3', 'y0', 'y1', 'y2', 'y3')
ndim = len(pnames)

# the interderometric data
infile = open(drk + 'dres.pkl','rb') #
results = pickle.load(infile)
infile.close()

weights = np.exp(results['logwt'] - results['logz'][-1])    # weighting the each nested 
postsamples = resample_equal(results.samples, weights)      # sample for posterior distribution after weighting the each nested 
                                                            
# the percentile parameters value 16%, 50%, 84% 
med = np.zeros(ndim)
for i in range(ndim):
    p = np.percentile(postsamples[:, i], [16, 50, 84])
    med[i] = p[1]
    print('%s in range %7.16e %7.16e %7.16e' % (pnames[i],p[0],p[1],p[2]))

# plotting with posterior samples 
tex = dict(ha='center', fontsize=12)
text = dict(ha='center', fontsize=12)
fig = corner.corner([postsamples[:,0], postsamples[:,1]], labels= pnames[:2], color = 'g', truths=truth[:2], plot_datapoints=False,
                   label_kwargs=tex, title_kwargs= text, show_titles=True, title_fmt=".2e")


plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.08, hspace=0.15)
plt.savefig(dfig + "parameters.png") 
plt.show()

