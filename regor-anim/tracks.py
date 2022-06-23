import numpy as np
import matplotlib.pyplot as pl

pl.style.use('dark_background')

from globi import flags

from layout import Layout

mirzam = Layout()
mirzam.set_radec((6,22,42),(-17,57,21))
mirzam.set_latlon((28,45,25.79),(-17,53,14.39)) # CTA-N + MAGIC                          

adhara = Layout()
adhara.set_radec((6,58,37.6),(-28,58,19))
adhara.set_latlon((28,45,25.79),(-17,53,14.39)) # CTA-N + MAGIC  

source = adhara

x = []
y = []
a = []

diam_mst = 12
diam_magic = 17
diam_lst = 23

if '_D25_' in flags:
    fil = open('D25.txt')

    
nfixed = len(a)-1
print(len(a),'telescopes')

jd0 = 2460000
jd = jd0 + np.arange(96)/96

ua = []
va = []
wa = []
aa = []

u,v,w = source.get_uvw(jd,0,0,1)
jd = jd[w>0]

for i in range(len(a)):
    for j in range(i):
        u,v,w = source.get_uvw(jd,x[i]-x[j],y[i]-y[j],0)
        oa = np.ones(len(jd)) * (a[i]*a[j])
        ua.append(u)
        va.append(v)
        wa.append(w)
        aa.append(oa)


    
track_u = np.concatenate(ua)
track_v = np.concatenate(va)
track_a = np.concatenate(aa)

if __name__ == '__main__':

    pl.xlim((-700,700))
    pl.ylim((-500,500))
    
    Amax = max(np.concatenate(aa))
    
    pl.title('Baseline in metres')
    pl.gca().set_aspect('equal')
    for b in range(len(aa)):
        #print(aa[b][0])
        if b >= nfixed*(nfixed-1)//2:
            col = 'magenta'
        else:
            col = 'white'
        pl.plot(-ua[b],-va[b],linewidth=aa[b][0]/50,color=col)
        pl.plot(ua[b],va[b],linewidth=aa[b][0]/50,color=col)
        '''
        if aa[b][0] == diam_lst**2:
            ls = 'dashed'
        elif aa[b][0] == diam_lst*diam_mst:
            ls = 'dashdot'
        elif aa[b][0] == diam_mst**2:
            ls = 'dotted'
        else:
            ls = 'solid'
        pl.plot(-ua[b],-va[b],linestyle=ls,color=col)
        pl.plot(ua[b],va[b],linestyle=ls,color=col)
        '''
    pl.tight_layout()
    #pl.savefig('baselines_'+flags.split('_S')[0]+'.png')
    pl.show()

