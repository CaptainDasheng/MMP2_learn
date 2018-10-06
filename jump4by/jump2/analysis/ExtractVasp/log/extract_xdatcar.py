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
#  learning.
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
Module for extract data from vasp output:XDATCAR. 
"""
import os
import numpy as np
from extract_basics import BaiscParam

__contributor__ = 'Xin He, Xingang Zhao'
__edited_date__ = '2017-04-29'

class get_ProStructure(BasicParam):
    """
    Class to get the structures of the trajectory from the XDATCAR 
            during the optimization/scf/nonscf by VASP calculation, 
            and append those structures into a lab.

    return:: Struct{'S1':{'structure':Structure, 'Energy':float}}
    """

    def __init__(self, path=None, *args, **kwargs):
	"""
	initialize the pathway
	"""
        
	super(get_ProStructure, self).__init__()

	# structure % 
        if path is None: path = os.getcwd()
        self.__path = path
   
    @property
    def structures(self):
	self.__read_structures()
        
	return self.__structure

    def __read_structures(self):

        """
        Read in the data in XDATCAR.
        """
	
	from jump2.structure import Structure

	structure = None
	from os.path import join 

	path = join(self.__path, 'XDATCAR')

	isif = self.isif  # 3 for relax cell, 2 for relax ions %  
	
	with open(path, 'r') as f: # open XDATCAR %
	    while True:	
	        line = f.readline()
                if not line: break
                if isif >= 3:
                    name = line # system name % 
		    scale, lattice, species, numbers =\
		    self.__read_lattice(f)

		if isif <3 and not structure:
                    name = line # system name % 
                    scale, lattice, species, numbers =\
                    self.__read_lattice(f)

		# structure of each step % 
                if line.startswith('Direct configuration='):
    
                    count = str('S'+line.split('=')[-1]).strip()
                    positions = self.__read_position(f, sum(numbers))
    	            
                    if not structure: structure = {}

                    structure[count] = Structure()
                    structure[count].comment   = name
                    structure[count].scale     = scale
                    structure[count].elements  = np.array(species)
                    structure[count].numbers   = np.array(species)
                    structure[count].lattice   = np.array(lattice)
                    structure[count].positions = np.array(positions)
                    structure[count].contraints= None
    
    def __read_lattice(self, file=None):
	
        line = file.readline()
        scale = float(line) # lattice scale % 
        
        # lattice % 
        for i in xrange(0,3):
            line = file.readline()
            lattice.append([float(v) for v in line.split()])
        
        # species in structure % 	
        line = file.readline()
        species = [str(s) for s in line.split()]
        
        # numbers of each element % 
        line = file.readline()
        numbers = [int(n) for n in line.split()]
    
	return scale, lattice, species, numbers
	
    def __read_position(self, file=None, num=None):
	"""
	return one structure position.
        """

        count = 0 
        position = []
        while count < num:
            line = file.readline()
            position.append([float(p) for p in line.split()])
            count += 1
	
	return position 
	#for i in xrange(0,5):
########name  = lines[0]
########scale = lines[1]
########for i in xrange
########  
########totalines = len(lines)
########
########elements = [str(s0) for s0 in lines[5].split()]
########numbers = [int(s1) for s1 in lines[6].split()]
########totalNumbers = sum(item for item in numbers)
########type = lines[14].split()[0]
#######
########if lines[7].split() == []:
########    count = 1+totalNumbers
########    step = (len(lines)-7)/(totalNumbers+1)
########    #lattice
########    lattice = []
########    for j in range(2,5):
########        temp = np.array([float(S0) for S0 in lines[j].split()])
########        if (lattice == []):
########            lattice = temp
########        else:
########            lattice = np.vstack([lattice,temp])
########    #positions
########    for i in range(int(step)):
########        positions = []
########        for k in range(8+i*count,8+i*count+totalNumbers):                   
########            temp = np.array([float(S0) for S0 in lines[k].split()])
########            if (positions == []):
########                positions = temp
########            else:
########                positions = np.vstack([positions,temp])                       


########        #dictionary
########        pdict = {}
########        pdict['comment'] = 'none'
########        pdict['elements'] = elements
########        pdict['positions'] = positions
########        pdict['lattice'] = lattice
########        pdict['numbers'] = numbers
########        pdict['type'] = type
########        pdict['constraints'] = 'none'
########        print pdict                    


########
########else:
########    count = 8+totalNumbers
########    step = lines[-totalNumbers-1].split()[-1]
########    for i in range(int(step)):
########        #lattice
########        lattice = []
########        for j in range(9+i*count,12+i*count):
########            temp = np.array([float(S0) for S0 in lines[j].split()])
########            if (lattice == []):
########                lattice = temp
########            else:
########                lattice = np.vstack([lattice,temp])

########                #position
########                positions = []
########        for k in range(15+i*count,15+i*count+totalNumbers):
########            temp = np.array([float(S0) for S0 in lines[k].split()])
########            if (positions == []):
########                positions = temp
########            else:
########                positions = np.vstack([positions,temp])
########
########        #dictionary
########        pdict = {}
########        pdict['comment'] = 'none'
########        pdict['elements'] = elements
########        pdict['positions'] = positions
########        pdict['lattice'] = lattice
########        pdict['numbers'] = numbers
########        pdict['type'] = 'none'
########        pdict['constraints'] = 'none'
########        print pdict 

                     
#xdatcar().read_xdatcar()
