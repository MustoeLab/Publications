##############################
#   AUROC and  ESI Plotter   #
#                            #
#     Anthony M. Mustoe      #
#          (c) 2023          #
##############################


import numpy as np
import matplotlib.pyplot as plot
from sklearn.metrics import roc_curve, roc_auc_score


def structureinf(r, p, u):
    
    si = 0
    for x in range(len(r)):
        if p[x] > 0:
            xp = p[x]/(p[x]+u[x])
            si += r[x]*xp*np.log2(xp)

        if u[x] > 0:
            xu = u[x]/(p[x]+u[x])
            si += r[x]*xu*np.log2(xu)
    
    return 1+si


def computeRatio(react, pairstatus, bins=20):

    ss, ed = np.histogram([react[i] for i,v in enumerate(pairstatus) if v==1], range=(0,1), bins=bins)
    ss = np.array(ss, dtype=float)/np.sum(ss)
    
    bp, ed = np.histogram([react[i] for i,v in enumerate(pairstatus) if v==0], range=(0,1), bins=bins)
    bp = np.array(bp, dtype=float)/np.sum(bp)
    
    r, ed = np.histogram(react, range=(0,1), bins=20)
    r = np.array(r, dtype=float)/np.sum(r)

    ratio = []
    for i in range(len(ss)):
        ratio.append(ss[i]/bp[i])
    
    si = structureinf(r,bp,ss)

    return ed[:-1], ratio, si




x = []
y = []
z = []
n = []
b = []

for i in range(100000):
    r = np.random.random()
    
    # single-stranded
    if r>0.5:
        b.append(1)
        if np.random.random()>0.5:
            x.append(np.random.beta(2,2))
        else:
            x.append(np.random.beta(6,4))
        
        if np.random.random()>0.5:
            y.append(np.random.beta(10,1))
        else:
            y.append(np.random.beta(1,5))
        
        z.append(np.random.beta(10,1))
        n.append(np.random.beta(1,5))

    # base-paired
    else:
        x.append(np.random.beta(1,2))
        y.append(np.random.beta(1,5))
        z.append(np.random.beta(1,10))
        n.append(np.random.beta(1,5))
        b.append(0)
    

fig,ax = plot.subplots(2,4)

ax[0,0].hist([x[i] for i,v in enumerate(b) if v==0], bins=50, histtype='step', color='C0', density=True)
ax[0,0].hist([x[i] for i,v in enumerate(b) if v==1], bins=50, histtype='step', color='C0', alpha=0.7, density=True)
ax[0,1].hist([y[i] for i,v in enumerate(b) if v==0], bins=50, histtype='step', color='C1', density=True)
ax[0,1].hist([y[i] for i,v in enumerate(b) if v==1], bins=50, histtype='step', color='C1', alpha=0.7, density=True)
ax[0,2].hist([z[i] for i,v in enumerate(b) if v==0], bins=50, histtype='step', color='C2', density=True)
ax[0,2].hist([z[i] for i,v in enumerate(b) if v==1], bins=50, histtype='step', color='C2', alpha=0.7, density=True)
ax[0,3].hist([n[i] for i,v in enumerate(b) if v==0], bins=50, histtype='step', color='C3', density=True)
ax[0,3].hist([n[i] for i,v in enumerate(b) if v==1], bins=50, histtype='step', color='C3', alpha=0.7, density=True)



ax[1,0].plot([0, 1], [0, 1], color='k', linestyle='--')


fpr, tpr, thr = roc_curve(b, x, drop_intermediate=False)
ax[1,0].plot(fpr, tpr, label='D1 AUC={:.2f}'.format(roc_auc_score(b,x)))
fpr, tpr, thr = roc_curve(b, y, drop_intermediate=False)
ax[1,0].plot(fpr, tpr, label='D2 AUC={:.2f}'.format(roc_auc_score(b,y)))
fpr, tpr, thr = roc_curve(b, z, drop_intermediate=False)
ax[1,0].plot(fpr, tpr, label='D3 AUC={:.2f}'.format(roc_auc_score(b,z)))

fpr, tpr, thr = roc_curve(b, n, drop_intermediate=False)
ax[1,0].plot(fpr, tpr, label='D4 AUC={:.2f}'.format(roc_auc_score(b,n)))

ax[1,0].legend()



ax[1,1].plot((0,1),(1,1), color='k', linestyle='--')

edge, ratio, si1 = computeRatio(x,b)
print('{} SI = {:.2f}'.format('D1', si1))
ax[1,1].plot(edge, ratio, label='D1')

edge, ratio, si2 = computeRatio(y,b)
print('{} SI = {:.2f}'.format('D2', si2))
ax[1,1].plot(edge, ratio, label='D2')

edge, ratio, si3 = computeRatio(z,b)
print('{} SI = {:.2f}'.format('D3', si3))
ax[1,1].plot(edge, ratio, label='D3')

edge, ratio, si4 = computeRatio(n,b)
print('{} SI = {:.2f}'.format('D4', si4))
ax[1,1].plot(edge, ratio, label='D4')

ax[1,1].set_yscale('log')


ax[1,2].bar((1,2,3,4), (si1,si2,si3,si4), color=('C0','C1','C2','C3'))
ax[1,3].bar((1,2,3,4), (roc_auc_score(b,x),roc_auc_score(b,y),roc_auc_score(b,z),roc_auc_score(b,n)), color=('C0','C1','C2','C3'))

outPath = 'Figure_S5_Base.pdf'
fig.savefig(outPath, dpi=100)


