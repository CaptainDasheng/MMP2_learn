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
        self.efermi = None
	
	self.__get_efermi__()	
	self.__get_eigenvalue__()

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
	self.__kpoints__ = kpoints 

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

	plotBS = True
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

    	        
    # claculate the distance of two kpoints % 
    def calcKpointsDistance(self, kpoints):
        
        kpointsDistance = []
        for i in range(0, kpoints.shape[0]):
            if (i == 0):
                kpointsDistance = np.array([0.0])
            else:
                c1 = kpoints[i]
                c0 = kpoints[i-1]

                a, b = self.getLatticeVector()
                dx = (c1[0]-c0[0])*b[0][0] + (c1[1]-c0[1])*b[1][0] + (c1[2]-c0[2])*b[2][0]
                dy = (c1[0]-c0[0])*b[0][1] + (c1[1]-c0[1])*b[1][1] + (c1[2]-c0[2])*b[2][1]
                dz = (c1[0]-c0[0])*b[0][2] + (c1[1]-c0[1])*b[1][2] + (c1[2]-c0[2])*b[2][2]
                delta = math.sqrt(dx*dx + dy*dy + dz*dz)

                temp = kpointsDistance[i-1] + delta
                kpointsDistance = np.hstack([kpointsDistance, temp])
        
        return kpointsDistance

    # get high symmetry points from KPOINTS
    def getHighSymmetryPoints(self):
        infile = open(self.path + "/KPOINTS")
        
        infile.readline() # comment line
        string = infile.readline()
        nhspoints = int(string)

        infile.readline() # skip
        string = infile.readline()

        # kpoints and kpointsName
        hspoints = []
        hspointsName = []
        if (string.startswith("reciprocal")):

            while(string):
                string = infile.readline()
                if (string != ""):
                    temp = np.array([float(s0) for s0 in string.split()[:3]])
                    name = string.split()[-1]
                    if (name == "\Gamma" or name == "\gamma"):
                        name = "$\Gamma$"
                    name = np.array([name])
                        
                    if (hspoints == []):
                        hspoints = temp
                        hspointsName = name
                    else:
                        if (hspointsName[-1] != name):
                            hspoints = np.vstack([hspoints, temp])
                            hspointsName = np.hstack([hspointsName, name])
        
        # high symmerty point distance
        kpoints, bands, xbar= self.getBand()
        distance = self.calcKpointsDistance(kpoints)

        hspointsDistance = np.array(distance[0]) # origin point
        for i in range(1, len(hspointsName)):
            index = nhspoints*i - 1
            hspointsDistance = np.hstack([hspointsDistance, distance[index]])

        return hspointsName, hspointsDistance, hspoints

    # plot band
    def plotBand(self):
        kpoints, bands, xbar = self.getBand()
        distances = self.calcKpointsDistance(kpoints)
        
        # plot bands
        for i in range(0, bands.shape[0]):
	    bands[i] =bands[i]+ xbar
            plt.plot(distances, bands[i], "b")

        # plot high symmerty lines
        hspointsName, hspointsDistance, hspoints = self.getHighSymmetryPoints()

        for i in range(0, len(hspointsDistance)):
            plt.axvline(x = hspointsDistance[i], color="k")

        plt.axhline(y = 0, linestyle="--") # fermi line

        plt.xticks(hspointsDistance, hspointsName)

        plt.ylim(-4,8)
        plt.xlim(distances[0], distances[-1])
        plt.tight_layout()
        plt.savefig("BandStructure.png", dpi=300)

