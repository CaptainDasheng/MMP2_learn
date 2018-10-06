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
Module for extract data from vasp output: data for band structure. 
"""

class BandStructure(BasicParam):

    """
    Get basic data associated with band structure

    tribute:
        vbm::tuple, (40,(a,b,c), -0.5)
        cbm::tuple, (41,(a,b,c),  0.5)
        ispin::integer, 1 or 2
        nelect:: integer, total number of electron; 
        symbols:: high symmetrized symbols;
        eigenvalue:: list, eigenvalues at each band;
        kpoints:: list, organized kpoints;
        gap:: dict, {'direct':{'x':(0.,0.0.),'value':1.0},
                     'indirect':{'x':(0.,0.,0.), 'value':0.5}}
        projected_bands:: dict of projected orbitals;
    """

    def __init__(self, path=None, symmetry_path=None, projected=False, *args, **kwargs):
        
        super(BandStructure, self).__init__()

        if path is None:
            path = os.path.getcwd()
        self.path = path 


    @property
    def vbm(self):
        """
        eigenvalue of valence band maximum.

        return::
            tuple(num,         # the num of points %    
                 [0., 0. 0.],  # coordination %  
                 -1.5)         # eigenvalue % 
        """
        return None 

    @property
    def cbm(self):
        """
        eigenvalue of conduction band minimum.

        return::
            tuple(num,         # the num of points %    
                 [0., 0. 0.],  # coordination %  
                 -1.5)         # eigenvalue % 
        """
        return None

    @property
    def ispin(self):
        """spin==1 or 2"""
        return None

    @property
    def efermi(self):
        #fermi_level % 
        return None

    @property
    def gap(self):
        # direct gap value % 
        return None

    @property
    def nelect(self):
        """the total number of electrons"""
        return None

    @property
    def symbols(self):
        # symbols of high symmetrized points % 
        return None

    @property 
    def eigenvalue(self):
        # eigenvalues of all the bands %
        return None 

    @property
    def kpoints(self):
        # kponts along high symmetrized pathway % 
        return None 

    @property
    def projected_bands(self):

        # projected bands % 
        return None

    def __iterm_path__(self):
        return None 
