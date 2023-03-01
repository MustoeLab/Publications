############################################
#        Plot G to A Mutation Rate         #
#                                          #
#           David Mitchell III             #
#                 (c) 2023                 #
############################################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys
import RNAtools2 as RNAtools
import numpy as np
from ReactivityProfile import ReactivityProfile


xtiks = (1,1.1,1.2)
labels = ('Unbuffered', 'pH 7.2', 'pH 8')

## Creating Plots ##
names = ('AUC', 'Mutation Rate')
fig, ax = plot.subplots(2,1, figsize=(3,4))
plot.rc('xtick', labelsize=10)
plot.rc('ytick', labelsize=10)

data_use = (.00355, .00329, .00274)
ylimrange = (0.0001, 0.1)
yscale_use = 'log'

ax[0].plot(xtiks, data_use, 'o', ms = 8 , mec = 'gray', mfc = 'None', ls = '-')
ax[0].set_ylim(ylimrange)
ax[0].set_ylabel('Mutation Rate')
ax[0].set_yscale(yscale_use)

plot.tight_layout()

## Creating Output PDF File ##
outPath1 = 'Output/RMRP_pH_Titration_GA_Mut_Rate.pdf'
plot.savefig(outPath1, dpi=100, bbox_inches="tight", transparent=True)

print('\nCompleted plot generation.\n')