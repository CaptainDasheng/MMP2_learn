#!/usr/bin/env python
#coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import os
import math

class BandStructure(object):

    """
    Get the eigenvalues from OUTCAR/EIGENVAL/vasprun.xml, default EIGENVAL files;

    method:
	1) band_gap: get direct/indirect band gap, for instance, 
			 'direct': VC: (0.,0.,0.), CC: (0., 0.,  0.), 1.0 eV;
		       'indirect': VC: (0.,0.,0.), CC: (0., 0.5, 0.), 1.0 eV;
	
	2) get_efermi: get fermi level;
	
	3) plot_bandStructure: plot band structure;

			 
    """
    def __init__(self, path=None, *args, **kwargs):

	if path is None: path = os.getcwd()
        self.__path__ =  path
	self.__nume__ = None
	self.__ispin__ = None 
        self.efermi = None

    # get fermi energy level from OUTCAR %
    def __get_efermi__(self):
	"""
	"""

        path = os.path.join(self.__path__, 'OUTCAR')
        string = os.popen("grep 'E-fermi' {0}".format(path)).readline()
        efermi = float(string.split()[2])
        
	self.efermi = (efermi, 'eV')

    # get band from EIGENVAL %
    # kpoints: n*3 array %
    # bands: nbands*nkpoints array %
    def __get_eigenvalue__(self):

        efermi = self.efermi[0]
        
	path = os.path.join(self.__path__, 'EIGENVAL')

	with open(path, 'r') as f:
	    lines = f.readlines()
       
	#ispin %  
	ispin = int(lines[0].split()[-1])
	
        # get total number of kpints and bands%
	num_e, num_k, num_b = [int(i) for i in lines[5].split()]

	self.__nume__ = num_e
	self.__ispin__ = ispin

	if ispin == 1:
	    num_e = int(round(num_e/2.0))
	
        # read bands
        bands = []
        kpoints = []
	
	direct_gap = 0.
	indirect_gap = 0.

	vbm = -1000.0
	cbm =  1000.0
		
	for k in xrange(0, num_k):
	    line = lines[6+k*(num_b+2):6+(k+1)*(num_b+2)]

	    each_k = tuple([ float(v) for v in line[1].split()[0:3]])
	    kpoints.append(each_k)

	    each_b = np.array([float(e.split()[-1]) for e in line[2:]])
	    bands.append(each_b - efermi)
	
	    if each_b[num_e] <= cbm: 
		cbm = each_b[num_e]
		cbm_k = each_k			    
	    if each_b[num_e -1] >= vbm: 
		vbm = each_b[num_e-1]
		vbm_k = each_k
                direct_gap = each_b[num_e] - each_b[num_e -1]

	indirect_gap = cbm - vbm 

	# gap % 
	direct = False
	if cbm_k == vbm_k: 
	    direct = True 
	transition = "VBM:({0})->CBM:({1})".format(\
		      ', '.join( str(l) for l in list(vbm_k)),\
		      ', '.join( str(l) for l in list(cbm_k)))

	self.__gap__ = {'is_direct':direct, \
			'transition':transition, \
			'direct_gap':(direct_gap, 'eV'),\
			'indirect_gap':(indirect_gap, 'eV')}
	
	self.__eigvalues__ = np.array(bands).T
	self.__kpoints__ = np.array(kpoints) 

    def get_bandstructure(self, *args, **kwargs):
        self.__get_efermi__()
        self.__get_eigenvalue__()

    @property
    def show_figure(self):

        if self.__figure__ is not None:
            self.__figure__.show()
        else:
            print "Run plotBandStructure first"

    @property
    def save_figure(self, name='bandstrucutre.png', **kwargs):
        if self.__figure__ is not None:
            self.__figure__.save(name, **kwargs)
        else:
            print "Run plotBandStructure first"
        
    @property 
    def ispin(self):
	return self.__ispin__

    @property
    def nelect(self):
	return self.__nume__

    @property
    def band_gap(self):
	"""
	detailed info for band gap.
	"""
	return self.__gap__

    @property
    def direct_gap(self):
	"""
	direct band gap and transition info.
	"""

	return self.__gap__['direct_gap'], self.__gap__['transition']

    @property 
    def indirect_gap(self):
	"""
	direct band gap and transition info.
	"""

	return self.__gap__['indirect_gap'], self.__gap__['transition']


    def __BS_input__(self):

	"""
	proccess the kpoints for plot band structure.
	"""	
	plotBS = None
	path = os.path.join(self.__path__, 'OUTCAR')
	
	with open(path, 'r') as f:
	    lines = f.readlines()
	
	for line in lines:
   	    r = re.compile(r'\s*k-points\s*in\s*reciprocal\s*lattice\s*and\s*weights:\s*')
            found = re.search(r, line)
            if found:
		if 'INCAR' in line or 'KPOINTS' in line:
		    plotBS = False
		elif 'symmetry lines' in line:
		    plotBS = True

	kpoints = self.__kpoints__ 

	xlabel = [(0, kpoints[0])]
	line_x = [0.0]
	# reorgnize the kpoints % 
	if plotBS is not None:
	    if plotBS is True:
		for i in xrange(0,len(kpoints)-1):
		   line_x.append(np.linalg.norm(np.array(kpoints[i+1])- np.array(kpoints[i])))
		   if kpoints[i] == kpoints[i+1]:
			xtlabel.append((i, kpoints[i]))
		xtlabel.append((i+1,kpoints[i+1]))
	
	    # 3D plot (under developed) % 
	    elif plotBS is False:
		kpoints = kpoints.T
		x,y  = [kpoints[0], kpoints[1]]
		[x_input, y_input] = np.meshgrid(x,y)	
		pass
	
	return (line_x, xlabel)

    # plot band
    def plot_bandstructure(self, ylim=-4.0, ymax=6.0, cbm_shift=0.0, 
				xlabels=None, color='k', linewidth=2.0, 
				filled=False, *args, **kwargs):
       
	
	x, xlabel = self.__BS_input__() 
	eigenvalues = self.__eigvalues__
	
	if self.ispin == 1:
	    num_vb = int(round(self.nelect/2.0))
	else:
	    num_vb = self.nelect
	
	count = 0 
        # plot bands
        for y in eigenvalues:
	    if count == num_vb: y = np.array(y) + shift
            plt.plot(x, y, '-', color=color, linewidth=linewith)
	    count += 1

        # plot high symmerty lines
        for i in xlabel:
            plt.axvline(x = x[i[0]], color="b")

        plt.axhline(y = self.efermi, linestyle="--") # fermi line

        plt.xticks(hspointsDistance, hspointsName)

        plt.ylim(ylim, ymax)
        plt.xlim(x[0], x[-1])
        plt.tight_layout()
	
        self.__figure__ = plt 
        #plt.savefig("BandStructure.png", dpi=300)

