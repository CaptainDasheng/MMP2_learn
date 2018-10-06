#!/usr/bin/env python 

from matplotlib import pyplot as plt 
import numpy as np 


x=np.linspace(-5,5, 100)
y=np.power(x,2)

fig= plt.figure()
ax=plt.subplot(111)

#ax.fill_between(x, y, [25]*len(x))
#ax.arrow(4, 4, 0.5, 0.5)
#ax.annotate('testest', xy=(5,5), xytext=(50,50), xycoords='data',arrowprops=dict(arrowstyle='->'))
ax.annotate('max', xy=(2,1), xytext=(3,1.5), arrowprops=dict(facecolor='k', arrowstyle='->'))
plt.savefig('fi.png',dpi=300)
