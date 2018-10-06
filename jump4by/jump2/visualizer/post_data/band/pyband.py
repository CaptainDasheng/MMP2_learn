#!/usr/bin/env python
import os, numpy, math 
import cPickle as pickle
import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import sys
from matplotlib.ticker import MultipleLocator

def plotBandStructure(ax, indat):
    hkpt=[]
    symp=[]
    # get the data
    (Symbols,Numbands, Nkpoints,dis, bande) = pickle.load(open(indat, 'rb'))
    #plt.scatter(x1, y1, s= z1, color='b', marker='.')
    for kk in xrange(Numbands):
        ax.plot(dis, bande[kk], 'k-')
    #plt.ylim(-10.0, 10.0)
    
    ax.set_ylim(-5.0, 6.0)
    ax.set_xlim(dis[0], dis[-1])
    ax.axhline(linewidth=0.5, color='r', linestyle='--')
    ax.axes.linewidth=2.0
    ax.set_xlabel("Band Structure", fontsize=12)
    ax.set_ylabel("Energy (eV)", fontsize=12)
    
    for kk in range(len(Symbols)):
        if Symbols[kk] != '':
            hkpt.append(kk)
            symp.append(Symbols[kk])
    for ii in range(1,len(hkpt)-1):
        ax.axvline(dis[hkpt[ii]])
    Symlocs=[]
    Symbols=[]
    Symlocs.append(dis[0])
    Symbols.append(symp[0])
    for ii in xrange(len(hkpt)):
	if ii%2 != 0:
	    Symlocs.append(dis[hkpt[ii]])
 	    Symbols.append(symp[ii])	
    ax.set_xticks(Symlocs)
    for i in xrange(len(Symbols)):
	if 'Gamma' in Symbols[i]:
	    Symbols[i] = 'G' 
    ax.set_xticklabels(Symbols, fontsize=10)

#======main program ==========
fig = plt.figure(figsize=(6,4))
ax1= plt.subplot(121)
ax2= plt.subplot(122)
if os.path.exists('AgIn.data'):
    plotBandStructure(ax1, 'AgIn.data')
    plotBandStructure(ax2, 'ZnSn.data')
plt.tight_layout()
ax1.set_title('$\mathregular{Cs_{4}[AgIn]_{2}Cl_{12}}$')
ax2.set_title('$\mathregular{Cs_{4}[Ag_{2}ZnSn]Cl_{12}}$')
bdf   = 'BandStructure'
showf = False 
if len(sys.argv) == 2 and ( str(sys.argv[1]) == 'T' or str(sys.argv[1]) == 't'): 
    showf = True
elif len(sys.argv) >= 2:
    bdf = str(sys.argv[1])
if showf : 
    plt.show()
else:
    plt.savefig(bdf+".png")
