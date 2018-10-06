# coding: utf-8
# Copyright (c) JUMP2 Development Team.
# Distributed under the terms of the JLU License.


#=================================================================
# This file is part of JUMP2.
#
# Copyright (C) 2017 Jilin University
#
#  Jump2 is a platform for high throughput calculation. It aims to 
#  make simple to organize and run large numbers of tasks on the 
#  superclusters and post-process the calculated results.
#  
#  Jump2 is a useful packages integrated the interfaces for ab initio 
#  programs, such as, VASP, Guassian, QE, Abinit and 
#  comprehensive workflows for automatically calculating by using 
#  simple parameters. Lots of methods to organize the structures 
#  for high throughput calculation are provided, such as alloy,
#  heterostructures, etc.The large number of data are appended in
#  the MySQL databases for further analysis by using machine 
#  learning (under developed).
#
#  Jump2 is free software. You can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published 
#  by the Free sofware Foundation, either version 3 of the License,
#  or (at your option) and later version.
# 
#  You should have recieved a copy of the GNU General Pulbic Lincense
#  along with Jump2. If not, see <https://www.gnu.org/licenses/>.
#=================================================================

"""
This module defines the classes relating to visualize the optical properties.
"""

from __future__ import absolute_import 

__contributor__ = 'Xin-Gang Zhao'


import sys, re, os
import numpy as np
import xml.etree.cElementTree as et
from numpy import * 

class Optical(object):

    """
    cls CalcOptics:: calculate the relative properties that 
        relates to the complex dielectric matrix with frequency,

        including, 

    method::

        1) get_absorption: absorption constant;
        2) get_jointDOS: joint density of states;
        3) get_refractive: calculate the refractive index;
        4) get_energyloss: energy loss coefficient;
        5) get_extinction: extinction coefficient;
        6) get_birefringence: refractive index that depends on the 
           polarization and propagation direction of light;
        7) get_photoconductivity: (need more information);
	8) get_slme: predict the spectroscopic limited maximum
	   efficiency.
        9) __get_dielectric_matrix__: get the real/imaginary epsilon parts,
           (just for VASP output right now). 
       
    args:
        path:  the orginal data, default is current path;
        input: the file with real/imaginary espilon constants,
	       default is OUTCAR, and EPS/vaspun.xml can also be
	       operated.	
    """

    def __init__(self, path=None, input='vasprun.xml', smear=None,
		          shift=0.0, sigma=0.01, *args, **kwargs):
	
        if path is None: 
            path = os.getcwd()

	self.__smear__ = smear 
	self.__sigma__ = sigma
	self.__datain__ = None
	self.__dielectric__ = None

	self.__path__ = path
	self.__input__ = input 
	self.__get_dielectric__()

    
    # method dielectric matrix %
    def __get_dielectric__(self):

        """
        :: buit-in method to get the dielectric matrix, including, 

           hw:: frequency (energy points)
           real_epsilon: real part of the dielectric constant;
           imag_epsilon: imaginary part of the dielectric constant;

           epsilon = complex(real_epsilon, imag_epsilon)

        return {'hw':[],'epsilon':[]}
        """

	
	input = self.__input__
        path = os.path.join(self.__path__, input)

        # from vasp output OUTCAR % 
	if input == 'OUTCAR':
	    self.__outcar_epsilon__(path)

        # from vasp output vasprun.xml % 
	elif input == 'vasprun.xml':
	    self.__vasprun_epsilon__(path)

        # from post-processing output EPS % 
	elif input == 'EPS':
	    self.__eps_epsilon__(path)
        
        # from dict % 
        # elif input == 
	else:
	    try:
		self.__outcar_epsilon__(os.path.join(self.__path__, 'OUTCAR'))
	    except:
	        raise IOError("{0} is not valid format file type.".format(input))
    
    # get the real/imaginary dielectric constants from vasprun.xml %
    def __eps_epsilon__(self, path):
	
        eps_imag = []
        eps_real = []
        hw = []
	
	with open(path, 'r') as f:
	    lines = f.readlines()
	for line in linesi[1:]:
	    hw.append(float(line.split()[0]))
	    eps_imag.append([float(line.split()[-1])]*3 + [0., 0., 0.])
	    eps_real.append([float(line.split()[1])]*3 + [0., 0., 0.])

	self.__transform__(hw, eps_real, eps_imag)
	
    # get the real/imaginary dielectric constants from vasprun.xml %
    def __vasprun_epsilon__(self, path):

        eps_imag = []
        eps_real = []
        hw = []

        tree = et.ElementTree(file=path)
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
                                           eps_imag.append(np.array(tmp.split()))
                                           hw.append(float(tmp.split()[0]))
                       if l.tag == 'real':
                          for x in l.getchildren():
                              for v in x.getchildren():
                                   if v.tag == 'set':
                                      for d in v.getchildren():
                                          tmp = d.text
                                          eps_real.append(np.array(tmp.split()))
        
	self.__transform__(hw, eps_real, eps_imag)
	
    # get the real/imaginary dielectric constants from OUTCAR %
    def __outcar_epsilon__(self, path):	
	
	hw = []
	real_epsilon = []
	imag_epsilon = []
	real = False
        imaginary = False

	# get the hw/real/imaginary epsilon %	
        with open(path, 'r') as file:
            lines = file.readlines()
        for i in xrange(0, len(lines)):
            regex1 = re.compile(r'\s*frequency\s*dependent\s*IMAGINARY\s*DIELECTRIC\s*FUNCTION\s*')
            found = regex1.match(lines[i])
	
            # imaginary part % 
            if found:
                imaginary = True
                i = i + 3
                while imaginary:
                    if re.match(r'^\n|\n+(?=\n)|\n$|(\s)\n',lines[i]):
                        imaginary = False
                        i = i - 1
                    else:
                        data = lines[i].split()
			hw.append(float(data[0]))
			imag_epsilon.append([float(j) for j in data[1:]])
                    i += 1
    
            regex2 = re.compile(r'\s*frequency\s*dependent\s*REAL\s*DIELECTRIC\s*FUNCTION\s*')
            found = regex2.match(lines[i])

	    # real part % 
            if found:
                realnary = True
                i += 3
                while realnary:
                    if re.match(r'^\n|\n+(?=\n)|\n$|(\s)\n', lines[i]):
                        realary = False
			break
                    else:
                        data = lines[i].split()
			real_epsilon.append([float(j) for j in data[1:]])
                    i += 1
            i += 1
	self.__transform__(hw, real_epsilon, imag_epsilon)

    
    def __transform__(self, hw, real_epsilon, imag_epsilon):
	"""
	To transform the dielectronic constants, with/without broaden
	imaginary part, default is None, and Gaussian/Lorenz type can be used.
	"""
	# transform %     
	epsilon=[]
	
	dE = hw[1] - hw[0]
	for i in xrange(len(hw)):
	    # gaussian broaden % 
	    if smear == 'gaussian':
		val = ((len(hw) - i)*dE/sigma)**2/2.0
		if val < 15.0: imag_epsilon[i] = dE*exp(-val)*\
			np.array(imag_epsilon[i])/sqrt(2.0*pi)/sigma

	    # lorenz broaden %
	    if smear == 'lorenz':
		pass 
  	
	    xx = complex(real_epsilon[i][0],imag_epsilon[i][0])
	    yy = complex(real_epsilon[i][1],imag_epsilon[i][1])
	    zz = complex(real_epsilon[i][2],imag_epsilon[i][2])
	    xy = complex(real_epsilon[i][3],imag_epsilon[i][3])
	    yz = complex(real_epsilon[i][4],imag_epsilon[i][4])
	    zx = complex(real_epsilon[i][5],imag_epsilon[i][5])
	    
	    # diagnolization %  
            eps = linalg.eig(mat([[xx,        xy,       conj(zx)], \
		                  [conj(xy),  yy,             yz], \
                                  [zx,        conj(yz),       zz]]))
            epsilon.append(np.array(eps))

	self.__dielectric__ = {'hw': np.array(hw), 'epsilon': epsilon} 


    # method absorption %
    def get_absorption(self):
    #    from numpy import *
        """
        :: optical absorption constant derived from following formula,

	alpha = sqrt(2*[-real_epsilon +/- |real_espilon + imag_epsilon*j|]) 
	
        return a dict named absorption {'hw':[...], 'alpha':[...]}
	(averaged along 3 directions).
        """
        

	alpha = []
	hw = self.__dielectric__['hw']
	epsilon = self.__dielectric__['epsilon']

	eps_v = abs(epsilon)-np.real(epsilon)
 
	alpha = sum(hw*71618.96076*sqrt(eps_v.T))/3.0	

	return {'hw':hw, 'absoprtion':alpha}

    # method refragement % 
    def get_refragement(self):

	"""
        :: optical refractive index derived from following formula,

             n = sqrt(0.5*[|real_espilon + imag_epsilon*j| + real_epsilon])

        return a dict {'hw':[...], 'refragement':[...]}

	"""
	hw = self.__dielectric__['hw']
	epsilon = self.__dielectric__['epsilon']
	eps_v = abs(epsilon)+np.real(epsilon)
	n = sum(sqrt(0.5*(eps_v.T)))/3.0	

	return {'hw':hw, 'refragement':n}

    #method extinction cofficient % 
    def get_extincton_coeff(self):

	"""
        :: extinction cofficient derived from following formula,

             ec = sqrt(0.5*[|real_espilon + imag_epsilon*j|-real_epsilon])

        return a list named ext = [...]
	"""
	hw = self.__dielectric__['hw']
	epsilon = self.__dielectric__['epsilon']
	eps_v = abs(epsilon)-np.real(epsilon)

	ec = sum(sqrt(0.5*(eps_v.T)))/3.0	
	return {'hw':hw, 'refragement':c}

    # method energy loss coefficient %
    def get_energyloss(self):

    #    from numpy import *
        """
        :: energy loss coefficient derived from following formula,
        
             elc = imag_espilon/|real_espilon + imag_epsilon*j|^2

        return a list named elc = [...]
        """
	hw = self.__dielectric__['hw']
	epsilon = self.__dielectric__['epsilon']
	eps_v = np.imag(epsilon)/abs(epsilon)**2
	l = sum(eps_v.T)/3.0	

	return {'hw':hw, 'refragement':l}

    # method birefrigence %
    def get_birefrigence(self):

        """
        :: optical birefigence derived from following formula,
        
             biR = xxxx

        return a list named biR = [...]
        """
       	# under developed % 
        pass


    # method photo conductivity %
    def get_photoconductivity(self):

        """
        :: photo conductivity derived from following formula,
        
             p_conduct = |-j*(omiga/4pi)*(epsilon -1)|, 
	
	where epsilon = real_epsilon + imag_epsilon*j

        return a list named p_conduct = [...]
        """
        
	hw = self.__dielectric__['hw']
	epsilon = self.__dielectric__['epsilon']
	
	eps_v = abs(complex(0,-1)*hw/(4*pi)*(epsilon-1))
	pt = sum(eps_v.T)/3.0	

	return {'hw':hw, 'conduct':pt}

    # method birefrigence %
    def get_birefrigence(self):

        """
        :: optical birefigence derived from following formula,
        
             biR = xxxx

        return a list named biR = [...]
        """
       	# under developed % 
        pass


    # method photo conductivity %
    def get_photoconductivity(self):

        """
        :: photo conductivity derived from following formula,
        
             p_conduct = |-j*(omiga/4pi)*(epsilon -1)|, 
	
	where epsilon = real_epsilon + imag_epsilon*j

        return a list named p_conduct = [...]
        """
        
        

    # method SLME % 
    def get_SLME(self, thickness=6.0, shift=None, *args, **kwargs):
	"""
	:: get the spectroscopic limited maximum efficiency, 
           please cite:
		Yu, L.; Zunger A.; Phys. Rev. Lett. 2012, 108(6), 68701.
	return dict {'thickness':[...], 'efficiency':[...]}
        """

    # method SLME % 
    def get_SLME(self, thickness=6.0, shift=None, *args, **kwargs):
	"""
	:: get the spectroscopic limited maximum efficiency, 
           please cite:
		Yu, L.; Zunger A.; Phys. Rev. Lett. 2012, 108(6), 68701.
	return dict {'thickness':[...], 'efficiency':[...]}

	"""
	
	self.__slme_input__()

	data_in = self.__datain__
	
	if shift is not None:
	    data_in[0] += shift 
	    Eg += shift 

	ev = 1.60217648740E-19
	h = 6.626068E-34
	c = 299792458
	Pin = 1000/ev
	
	thick = []
	slme_v = []
	
	for L in np.linspace(0., thickness, thickness*100):
            Vc = 0.025851997434
            Isc = 0.0
            I0  = 0.0
            for l in xrange(len(data_in)-1) :
                hv0 = data_in[l][0]
                hv1 = data_in[l+1][0]
		#
                des = hv1 - hv0
                #
                aE0 = 1.0 - exp(-2.0*L*data_in[l][2])
                aE1 = 1.0 - exp(-2.0*L*data_in[l+1][2])

                is0 = data_in[l][1]*aE0
                is1 = data_in[l+1][1]*aE1

                Isc = Isc + (is0 + is1)*des/2.0
                
                irb0 = hv0**2/(exp(hv0/Vc)-1)*aE0
                irb1 = hv1**2/(exp(hv1/Vc)-1)*aE1
            #   irb0 = hv0**2/(exp(hv0/Vc)-1)*aE0*data_in[l][3]
            #   irb1 = hv1**2/(exp(hv1/Vc)-1)*aE1*data_in[l+1][3]
               
                I0 = I0 + (irb0 + irb1)*des/2.0

            I0 = I0 * 2.0*pi*c/(h*c/ev)**3 * exp(dEg/Vc)

            # calculate max IV = [ Isc - I0*e^{(V+dEg)/Vc} ] *V
            # dE = 0.001
            #npts = int(thres_Eg/0.001) 
            npts = int(Eg/0.001) 

            maxIV = 0
            for ll in xrange(npts) :
                Vap = ll*0.001
                IVtmp = Vap*( Isc - I0*exp(Vap/Vc)) 
                if IVtmp > maxIV :
                   maxIV = IVtmp
                   #%Vm = Vap

            slme = maxIV/Pin*100.0
                
            slme_v.append(slme)
	
	return {'thick':thick, 'slme':slme_v}
	
    # prepare input for calculating the SLME %
    def __slme_input__(self):
       

        from stand_am15 import am15 
 	
	alpha = self.get_absorption()
	refra = self.get_refragement()
	
	data = [alpha['hw'], alpha['absorption'], refra['refragment']**2] 

	data_in = []

	for l in range(1, len(am15)) :
	    x = am15[l].split()
	    hv  = float(x[0])
	    nhv = float(x[1])

	    for ll in range(len(alpha)-1) :
		if alpha[ll][0] <= hv and alpha[ll+1][0] >= hv :
		   fact = (hv - alpha[ll][0])/(alpha[ll+1][0] - alpha[ll][0])
		   tmp1 = alpha[ll][1]*(1-fact) + fact*alpha[ll+1][1]
		   tmp2 = alpha[ll][2]*(1-fact) + fact*alpha[ll+1][2]
		   data_in.append([hv, nhv, tmp1, tmp2])
		   break

	self.__datain__ = data_in
	del data_in, alpha, refra 

    def get_TaucPlot(self, r=1/2, *args, **kwargs):
	
	"""
	args r:
          r = 1/2,     direct allowed transition;
          r = 2,     indirect allowed transition;
          r = 3/2,   direct forbidden transition;
          r = 3,   indirect forbidden transition;
       
        See detailed derivations in:
        1) Optical Properties and Electonic Structure of Amorphous Germanium,
	   J. Tauc, R. Grigorovci and A. Vancu, Phys. Stat. Sol. 15, 627 (1966).
	2) Electronic Structure Optical Properties of Amorphous Ge and Si,
           Y. F. Tsay and D. K. Paul, Phys. Rev. B, 6, 2827 (1973).

	"""

	data = get_absorption()
	
	hw = data['hw']
	alpha = data['absorption']**r

	#fit a line (under developed) % 
	line = None	
	return {'hw': hw, 'tauc':alpha, 'line':line}


    # momentum matrix % 
    def get_momentum(self):
        
        """
        :: get the transition matrix derived from following formula,
        """
	
       	energy = []
	amp = []
	 
	path = os.path.join(self.__path__, 'MATRIX')
        if os.path.exists(path):
            pass
        else:
            raise AttributeError("File MATRIX doesn't exist, please check it out ")

	with open(path, 'r') as file:
	    lines = file.readlines()
	
	for line in lines[3:]:
	    if line.startswith('E='):
		energy.append(float(lin.split()[2]))	
		amp.append(float(lin.split()[5]))	
	
	return {'hw': np.array(energy), 'amplitude':np.array(amp)}

    # method jonit density of states %
    def get_jointDOS(self):

        """
        :: joint density states (JDOS) derived from following formula,
        
             jdos = xxxx

        return a list named jdos = [...]
        """
	path = os.path.join(self.__path__, 'JDOS')
        if os.path.exists(path):
            pass
        else:
            raise AttributeError("File JDOS doesn't exist, please check it out ")

        with open(path, 'r') as file:
            lines = file.readlines()

        hw = []
        jdos = []
 
	for line in lines[1:]:
	    hw.append(float(line.split()[0]))
            jdos.append(float(line.split()[1]))

        return {'hw':np.array(hw), 'jdos':jdos}

