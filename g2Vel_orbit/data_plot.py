import numpy as np
import matplotlib.pyplot as plt
from obslen import base

# observation Julian day and signals
jd = np.load("data/jd.npy")
mdata = np.load("data/sig.npy")
sdata = np.load("data/sig_noise.npy")

# three dimensional baseline 
east = np.real(base)
north = np.imag(base)

#signal with noise
plt.rcParams.update({'font.size': 30})#,'font.weight': 'bold'})
fig, ax = plt.subplots(6,sharex=True)
#fig.subtitle('Signal with Julian Day')
ax[0].plot(jd, sdata[0], color='royalblue', label = (east[0], north[0]))
ax[0].legend(bbox_to_anchor=(0.85, 0.85))
ax[1].plot(jd, sdata[1], color='green', label =(east[1], north[1]))
ax[1].legend(bbox_to_anchor=(0.85, 0.85))
ax[2].plot(jd, sdata[2], color='red', label = (east[2], north[2]))
ax[2].legend(bbox_to_anchor=(0.85, 0.85))
ax[3].plot(jd, sdata[3], color='purple', label = (east[3], north[3]))
ax[3].legend(bbox_to_anchor=(0.85, 0.85))
ax[4].plot(jd, sdata[4], color='brown', label = (east[4], north[4]))
ax[4].legend(bbox_to_anchor=(0.85, 0.85))
ax[5].plot(jd, sdata[5], color='darkorange', label = (east[5], north[5]))
ax[5].legend(bbox_to_anchor=(0.85, 0.85))
fig.text(0.5, 0.04, 'Julian Observation Day', ha='center', fontsize=30)
fig.text(0.01, 0.5, 'Correlated Signal According to baseline with noise', va='center', rotation='vertical', fontsize=30)
plt.subplots_adjust(left=0.1, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.gcf().set_size_inches(18, 18)
plt.savefig("data/fig/nsignal.png")
plt.show()

# model data from each baseline
plt.rcParams.update({'font.size': 30})#,'font.weight': 'bold'})
fig, ax = plt.subplots(6,sharex=True)
#fig.subtitle('Signal with Julian Day')
ax[0].plot(jd, mdata[0], color='royalblue', label = (east[0], north[0]))
ax[0].legend(bbox_to_anchor=(0.85, 0.85))
ax[1].plot(jd, mdata[1], color='green', label =(east[1], north[1]))
ax[1].legend(bbox_to_anchor=(0.85, 0.85))
ax[2].plot(jd, mdata[2], color='red', label = (east[2], north[2]))
ax[2].legend(bbox_to_anchor=(0.85, 0.85))
ax[3].plot(jd, mdata[3], color='purple', label = (east[3], north[3]))
ax[3].legend(bbox_to_anchor=(0.85, 0.85))
ax[4].plot(jd, mdata[4], color='brown', label = (east[4], north[4]))
ax[4].legend(bbox_to_anchor=(0.85, 0.85))
ax[5].plot(jd, mdata[5], color='darkorange', label = (east[5], north[5]))
ax[5].legend(bbox_to_anchor=(0.85, 0.85))
fig.text(0.5, 0.04, 'Julian Observation Day', ha='center', fontsize=30)
fig.text(0.01, 0.5, 'Correlated Signal According to baseline', va='center', rotation='vertical', fontsize=30)
plt.subplots_adjust(left=0.1, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.gcf().set_size_inches(18, 18)
plt.savefig("data/fig/signal.png")
plt.show()

