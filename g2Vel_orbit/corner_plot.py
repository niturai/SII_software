import numpy as np
from dynesty.utils import resample_equal
import pickle
import matplotlib.pyplot as plt
import corner

# the radius of stars 
R_a = 39.4503                      # Radius of source A in light seconds
R_b = 13.9236                       # Radius of source B in light seconds

# position of orbit in sky with time
dAB = np.load("data/dAB.npy")
X = dAB[0,0]                  # X-position of orbit
Y = dAB[1,0]                  # Y-position of orbit

# parameters, which are estimated
truth = [R_a, R_b, X, Y]
print('parameters',truth)


#number of parameter to estimate
pnames = ('$R_a (ls) $', '$R_b (ls) $', 'X', 'Y')
ndim = len(pnames)

# the interderometric data
infile = open('data/dres.pkl','rb')
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
plt.savefig("data/fig/corner.png")
plt.show()

