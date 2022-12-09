import numpy as np
from dynesty.utils import resample_equal
import pickle
import matplotlib.pyplot as plt
import corner
import math
from intfer_rl_ld import R_a, R_b, R1, X, Y, l #  

# parameters, which are estimated
truth = [R_a, R_b, R1, X, Y, l] #
print('parameters',truth)


#number of parameter to estimate
pnames = ('$R_a (ls) $', '$R_b (ls) $', 'R1', 'X', 'Y', 'l') #
ndim = len(pnames)

# the interderometric data
infile = open('data/dres_rl_ld.pkl','rb') #
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
tex = dict(ha='center', fontsize=14,fontweight='bold', color='black')
text = dict(ha='center', fontsize=12,fontweight='bold', color='black')
fig = corner.corner(postsamples, labels= pnames, truths=truth, plot_datapoints=False,
                    label_kwargs= tex, title_kwargs= text, show_titles=True, title_fmt=".2e")

plt.gcf().set_size_inches(8, 8)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.11)
plt.savefig("data/fig/corner_rl_ld.png") #
plt.show()

