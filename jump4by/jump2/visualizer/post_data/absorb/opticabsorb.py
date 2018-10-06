#!/usr/bin/env python 
from numpy import *
from scipy import linalg
import sys
import re
import cPickle as pickle 
import numpy as np
import xml.etree.cElementTree as et
from matplotlib import pyplot as plt
import numpy as np
#plt.figure(1)

def plotoptics(ax,energy=None,showl=True, col=['b'],
		shift=0.0,tensor='Averaged Absorption',absor=None):

        energy = np.array(energy) + float(shift)
	if showl:
        	ax.semilogy(np.array(energy),np.array\
        	(absor),'-', color='b',linewidth=2.0, label=tensor)
	else:
        	plt.plot(np.array(energy),np.array\
        	(absor),'-', color='b',linewidth=2.0, label=tensor)

#=====program main ======
import os 
#plot the absorption
fig=plt.figure(figsize=(6,3))
ax1= plt.subplot(121)
ax2= plt.subplot(122)
inf=['AgIn.abs', 'ZnSn.abs']
axf=[ax1, ax2]
for i in xrange(0,2):
    (Energyp, Absorpv, Reflexd) = pickle.load(open(inf[i], 'rb'))
    plotoptics(axf[i], energy=Energyp,absor=Absorpv)
    ax = axf[i]
    ax.set_xlim(0,3.5)
    ax.set_ylim(100,1000000)
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Absorption (cm^-1)')
    #ax.legend(loc=2)
showf = False
showl = False
bdf   = 'optics'

if len(sys.argv) == 2 and ( str(sys.argv[1]) == 'T' or str(sys.argv[1]) == 't'):
    showf = True
elif len(sys.argv) == 2 and (str(sys.argv[1]) == 'L' or str(sys.argv[1]) == 'l'):
    #bdf = str(sys.argv[1])
    showl = True
    showf = True
ax1.set_title('$\mathregular{Cs_{4}[AgIn]_{2}Cl_{12}}$')
ax2.set_title('$\mathregular{Cs_{4}[Ag_{2}ZnSn]Cl_{12}}$')
#plt.show()
plt.tight_layout()
plt.savefig(bdf+'.png',dpi=300)
#main script.
