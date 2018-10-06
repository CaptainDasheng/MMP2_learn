#!/usr/bin/env python 
from numpy import *
from scipy import linalg
import sys
import re
import numpy as np
import xml.etree.cElementTree as et

class Absorb():
    
	def __init__(path=None):

		if path == None: self.path = os.getcwd()
	
	def  GetEpsilon(self):

		Epsimag = []
		Epsreal = []
		Energyp = []

		tree = et.ElementTree(file='vasprun.xml')
		root = tree.getroot()

		for m in root:
		   for c in  m.getchildren():
		       if c.tag == 'dielectricfunction':
		           for l in c.getchildren():
		               if l.tag == 'imag':
		                   for x in l.getchildren():
	        	               for v in x.getchildren():
		                           if v.tag == 'set':
	                	               for d in v.getchildren():
	       		        	    	   tmp = d.text
	                                	   Epsimag.append(np.array(tmp.split()))
	       			   	           Energyp.append(float(tmp.split()[0]))
		               if l.tag == 'real':
		                  for x in l.getchildren():
	        	              for v in x.getchildren():
	               		           if v.tag == 'set':
	                        	      for d in v.getchildren():
	       	            			  tmp = d.text
	                                	  Epsreal.append(np.array(tmp.split()))
		return Energyp, Epsreal, Epsimag
	
	def CalAbsorb(self):
            Energyp = []
	    Epsreal = []
	    Absorpv = []
	    Epsimag = [] 
	    Reflexd = []
	    Energyp, Epsreal, Epsimag = self.GetEpsilon()
	    for i in range(len(Epsimag)):
	       Energy=0.
	       ImagXX=0.
	       ImagYY=0.
	       ImagZZ=0.
	       RealXX=0.
	       RealYY=0.
	       RealZZ=0.
	       ImagXY=0.
	       ImagYZ=0.
	       ImagZX=0.
	       RealXY=0.
	       RealYZ=0.
	       RealZX=0.
		
	       Energy=float(Energyp[i])
	       ImagXX=float(Epsimag[i][1])
	       ImagYY=float(Epsimag[i][2])
	       ImagZZ=float(Epsimag[i][3])
	       RealXX=float(Epsreal[i][1])
	       RealYY=float(Epsreal[i][2])
	       RealZZ=float(Epsreal[i][3])
	       ImagXY=float(Epsimag[i][4])	
	       ImagYZ=float(Epsimag[i][5])	
	       ImagZX=float(Epsimag[i][6])	
	       RealXY=float(Epsreal[i][4])	
	       RealYZ=float(Epsreal[i][5])	
	       RealZX=float(Epsreal[i][6])	
	       
	       Cxx = complex(RealXX, ImagXX)
	       Cyy = complex(RealYY, ImagYY)
	       Czz = complex(RealZZ, ImagZZ)
	       Cxy = complex(RealXY, ImagXY)
	       Cyz = complex(RealYZ, ImagYZ)
	       Czx = complex(RealZX, ImagZX)
	       
	       C_eps = mat([[Cxx, Cxy, conj(Czx)],[conj(Cxy), Cyy, Cyz], [Czx, conj(Cyz), Czz]])
	       
	       eps_eig, eps_v = linalg.eig(C_eps)
	       hv = Energy
	       alpha_a1 = hv * 71618.96076 * sqrt(abs(eps_eig[0])-real(eps_eig[0]))
	       alpha_a2 = hv * 71618.96076 * sqrt(abs(eps_eig[1])-real(eps_eig[1]))
	       alpha_a3 = hv * 71618.96076 * sqrt(abs(eps_eig[2])-real(eps_eig[2]))
	       alpha_av = (alpha_a1 + alpha_a2 + alpha_a3)/3
	       Absorpv.append(alpha_av)
	
	       n1 = sqrt(0.5*(abs(eps_eig[0])+real(eps_eig[0])))
	       n2 = sqrt(0.5*(abs(eps_eig[1])+real(eps_eig[1])))
	       n3 = sqrt(0.5*(abs(eps_eig[2])+real(eps_eig[2])))
	       n_av = (n1+n2+n3)/3.0
	      
	       Reflexd.append(n_av)
	    return Energyp, Absorpv, Reflexd
	
	def savedata(self):
		import sys, os
		
		Energyp, Absorpv, Reflexd = self.CalAbsorb()
		sys.stdout = open('optics.data','w')
                print '%10s  %10s' % ('energy(ev)','absorb(cm^-1)')
                print
                for n in xrange(len(Energyp)):
                        print '%10.4f %10.4f' % (float(Energyp[n]),float(Absorpv[n]))
                sys.stdout.close()
	
	
	def plotoptics(self,energy=None,col=['b'],shift=0.0,tensor='Averaged Absorption',absor=None):
	        from matplotlib import pyplot as plt
	        import numpy as np
	        #plt.figure(1)
		showf = False
		showl = False
		bdf   = 'optics'
		if len(sys.argv) == 2 and ( str(sys.argv[1]) == 'T' or str(sys.argv[1]) == 't'):
		    showf = True
		elif len(sys.argv) == 2 and (str(sys.argv[1]) == 'L' or str(sys.argv[1]) == 'l'):
		    #bdf = str(sys.argv[1])
		    showl = True
		    showf = True

	        energy = np.array(energy) + float(shift)
		if showl:
	        	plt.semilogy(np.array(energy),np.array\
	        	(absor),'-', color='b',linewidth=2.0, label=tensor)
		else:
	        	plt.plot(np.array(energy),np.array\
	        	(absor),'-', color='b',linewidth=2.0, label=tensor)
	        plt.xlim(0,3.5)
#		plt.ylim(100,1000000)
	        plt.xlabel('Energy (eV)')
	        plt.ylabel('Absorption (cm^-1)')
	        plt.legend(loc=2)
	        plt.show()
	        plt.savefig(bdf+'.png',dpi=600)
#main script.
import os 
a=Absorb()
if not os.path.exists('optics.data'):
	a.savedata()

#plot the absorption
Energyp, Absorpv, Reflexd = a.CalAbsorb()
a.plotoptics(energy=Energyp,absor=Absorpv)
