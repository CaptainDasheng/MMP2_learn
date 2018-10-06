#!/usr/bin/env python
import os, numpy, math 
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np 

from jdos import *

fig=plt.figure(figsize=(8,4))
ax1=plt.subplot(121)
ax2=plt.subplot(122)
#plot JDOS % 
gap=[1.03,0.81]
mt=['AgIn.jdos', 'ZnSn.jdos']
fg=[ax1, ax2]
for i in xrange(0,2):
	ax = fg[i]
        PlotJDOS(ax, mt[i], shift=0.0)
        ax.set_ylim(0,6)
        ax.set_xlim(0,3)
        if i ==0: 
            ax.set_ylabel('JDOS (states/eV)\n', fontsize=12).set_weight('bold')
        ax.set_yticks([1,2,3,4,5,6])
        ax.set_yticklabels([1.0,2.0,3.0,4.0,5.0,6.0])
        ax.arrow(gap[i],0.5, 0.0, -0.2, fc='r', ec='r',
        	head_width=0.05, head_length=0.2)
        ax.text(gap[i]-0.5, 0.6, "$\mathregular{E_{g}="+str(gap[i])+"}$")
	ax.set_xlabel("$\mathregular{Energy \ (eV)}$")	
#plt.show()
ax1.set_title('$\mathregular{Cs_{4}[AgIn]_{2}Cl_{12}}$')
ax2.set_title('$\mathregular{Cs_{4}[Ag_{2}ZnSn]Cl_{12}}$')
plt.savefig('JDOS.png', dpi=300)
