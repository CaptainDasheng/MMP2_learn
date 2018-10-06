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


__author__ = ""
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import math
from matplotlib.font_manager import FontProperties

class BandStructure():

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
        
        #super(BandStructure, self).__init__()

        if path is None:
            path = os.getcwd()
        self.path = path

    @property
    def vbm(self):
        """
        eigenvalue of valence band maximum.

        return::
            tuple(num,         # the num of points %    
                 [0., 0., 0.], # coordination %
                 -1.5)         # eigenvalue % 
        """
        VBMtop = -100
        kpoint = 0
        bands = self.eigenvalue
        for i in range(0, bands.shape[0]):
            for j in range(0, bands.shape[1]):
                if bands[i][j] <= 0.5:
                    if bands[i][j] > VBMtop:
                        VBMtop = bands[i][j]
                        kpoint = j

        coordination = self.kpoints[kpoint]

        return (kpoint, coordination, VBMtop)

    @property
    def cbm(self):
        """
        eigenvalue of conduction band minimum.

        return::
            tuple(num,         # the num of points %    
                 [0., 0. 0.],  # coordination %  
                 -1.5)         # eigenvalue % 
        """
        CBMbottom = 100
        kpoint = 0
        bands = self.eigenvalue
        for i in range(0, bands.shape[0]):
            for j in range(0, bands.shape[1]):
                if bands[i][j] > 0.5:
                    if bands[i][j] < CBMbottom:
                        CBMbottom = bands[i][j]
                        kpoint = j

        coordination = self.kpoints[kpoint]

        return (kpoint, coordination, CBMbottom)

    @property
    def ispin(self):
        """spin==1 or 2"""
        temp = os.popen('grep ISPIN {0}'.format('OUTCAR')).readline()
        ispin = int(temp.split()[2])
        return ispin

    @property
    def efermi(self):
        #fermi_level %
        temp = os.popen("grep 'E-fermi' OUTCAR").readline()
        efermi = float(temp.split()[2])
        return efermi

    @property
    def gap(self):
        # direct gap value %

        direct_gap = 0.0
        direct_kpoint = []
        indirect_gap = 0.0

        cbm = self.cbm
        vbm = self.vbm

        if cbm[2] - vbm[2] > 10E-2:
            bands = self.eigenvalue
            direct_gap = 100.0
            for i in xrange(0, bands.shape[1]):
                vmaxer = -100.0
                cminer = 100.0
                for j in xrange(1, bands.shape[0]):
                    temp = bands[j][i]
                    if temp >= 0.0:
                        if temp < cminer:
                            cminer = temp
                    else:
                        if temp > vmaxer:
                            vmaxer = temp
                gap = cminer - vmaxer

                if gap < direct_gap:
                    direct_gap = gap
                    direct_kpoint = i

            indirect_gap = cbm[2] - vbm[2]
            if (cbm[1] == vbm[1]).all():
                isdirect =True

        return {'direct': {'x': self.kpoints[direct_kpoint], 'value': direct_gap},
                'indirect': {'x': self.kpoints[self.cbm[0]], 'value': indirect_gap}}

    @property
    def nelect(self):
        """the total number of electrons"""
        temp = os.popen("grep 'NELECT' OUTCAR").readline()
        nelect = int(float(temp.split()[2]))
        return nelect

    @property
    def symbols(self):
        # symbols of high symmetrized points %
        infile = open(self.path + "/KPOINTS")

        infile.readline()  # comment line
        string = infile.readline()
        nhspoints = int(string)

        infile.readline()  # skip
        string = infile.readline()

        # kpoints and kpointsName
        hspoints = []
        hspointsName = []
        if (string.startswith("reciprocal")):

            while (string):
                string = infile.readline()
                if (string != ""):
                    temp = np.array([float(s0) for s0 in string.split()[:3]])
                    name = string.split()[-1]
                    if (name == "Gamma" or name == "gamma"):
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
        kpoints = self.kpoints

        infile2 = open(self.path + "/POSCAR")
        infile2.readline()
        string = infile2.readline()
        latticeConstant = float(string)

        a = []  # lattice vector
        for i in range(0, 3):
            string = infile2.readline()
            temp = np.array([float(s0) * latticeConstant for s0 in string.split()])
            if (a == []):
                a = temp
            else:
                a = np.vstack([a, temp])

        volume = a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] + \
                 a[0][2] * a[1][0] * a[2][1] - a[0][0] * a[1][2] * a[2][1] - \
                 a[0][1] * a[1][0] * a[2][2] - a[0][2] * a[1][1] * a[2][0]

        b = np.zeros((3, 3))
        for i in range(0, 3):
            if (i == 0):
                j = 1
                k = 2
            elif (i == 1):
                j = 2
                k = 0
            else:
                j = 0
                k = 1
            c = np.zeros(3)
            c[0] = a[j][1] * a[k][2] - a[j][2] * a[k][1]
            c[1] = a[j][2] * a[k][0] - a[j][0] * a[k][2]
            c[2] = a[j][0] * a[k][1] - a[j][1] * a[k][0]
            for j in range(0, 3):
                b[i][j] = 2 * math.pi * c[j] / volume

        distance = []
        for i in range(0, kpoints.shape[0]):
            if (i == 0):
                distance = np.array([0.0])
            else:
                c1 = kpoints[i]
                c0 = kpoints[i - 1]

                dx = (c1[0] - c0[0]) * b[0][0] + (c1[1] - c0[1]) * b[1][0] + (c1[2] - c0[2]) * b[2][0]
                dy = (c1[0] - c0[0]) * b[0][1] + (c1[1] - c0[1]) * b[1][1] + (c1[2] - c0[2]) * b[2][1]
                dz = (c1[0] - c0[0]) * b[0][2] + (c1[1] - c0[1]) * b[1][2] + (c1[2] - c0[2]) * b[2][2]
                delta = math.sqrt(dx * dx + dy * dy + dz * dz)

                temp = distance[i - 1] + delta
                distance = np.hstack([distance, temp])



        hspointsDistance = np.array(distance[0])  # origin point
        for i in range(1, len(hspointsName)):
            index = nhspoints * i - 1
            hspointsDistance = np.hstack([hspointsDistance, distance[index]])

        return hspointsName, hspointsDistance, hspoints

    @property 
    def eigenvalue(self):
        # eigenvalues of all the bands %

        Efermi = self.efermi
        infile = open(self.path + "/EIGENVAL")

        # comment lines
        for i in range(0, 5):
            infile.readline()
        # get total number of kpints and bands
        string = infile.readline()
        totalKpoints = int(string.split()[1])
        totalBands = int(string.split()[2])

        # read bands
        bands = []
        counter = 0
        for i in range(0, totalKpoints):
            infile.readline()  # blank line
            infile.readline()  # kpoint line

            # get band value
            bk = []  # bands of a kpoint
            for j in range(0, totalBands):
                string = infile.readline()
                temp = np.array(float(string.split()[1]) - Efermi)
                if (bk == []):
                    bk = temp
                else:
                    bk = np.vstack([bk, temp])

            # bk->bands
            if (bands == []):
                bands = bk
            else:
                bands = np.hstack([bands, bk])
            counter += 1

        if (counter != totalKpoints):
            print "error! read number of kpoints", counter

        return bands

    @property
    def kpoints(self):
        # kponts along high symmetrized pathway %

        infile = open(self.path + "/EIGENVAL")

        # comment lines
        for i in range(0, 5):
            infile.readline()
        # get total number of kpints and bands
        string = infile.readline()
        totalKpoints = int(string.split()[1])
        totalBands = int(string.split()[2])

        # read bands
        kpoints = []
        for i in range(0, totalKpoints):
            infile.readline()  # blank line

            # get kpoint coordinate
            string = infile.readline()
            kp = np.array([float(s0) for s0 in string.split()[:3]])
            if (kpoints == []):
                kpoints = kp
            else:
                kpoints = np.vstack([kpoints, kp])

            for j in range(0, totalBands):
                infile.readline()

        return kpoints

    @property
    def projected_bands(self):
        # projected bands %

        infile = open(self.path + "/PROCAR")
        # comment lines
        infile.readline()
        # get total number of kpints and bands
        string = infile.readline()
        totalKpoints = int(string.split()[3])
        totalBands = int(string.split()[7])
        natom = int(string.split()[11])
        norbit = len(os.popen("grep -w ion PROCAR").readline().split()) -2

        # read bands
        ions = np.zeros((totalBands, totalKpoints, natom, norbit+1))
        counter = 0
        for i in range(0, totalKpoints):  # kpoint
            # skip blank line
            infile.readline()
            # skip kpoint line
            infile.readline()
            # get band value
            for j in range(0, totalBands):  # band
                # skip blank line
                infile.readline()
                # skip band line
                infile.readline()
                # skip blank line
                infile.readline()
                # skip comment line
                infile.readline()
                for k in range(0, natom):  # atom
                    string = infile.readline()
                    temp = np.array([s0 for s0 in string.split()])
                    for l in xrange(0, norbit+1):
                        ions[j][i][k][l] = float(temp[l + 1])
                # skip "tot..."
                infile.readline()
                if (j == totalBands - 1):
                    # skip blank line, enter next kpoint
                    infile.readline()

            counter += 1

        if (counter != totalKpoints):
            print "error! read number of kpoints", counter

        return ions

    def __iterm_path__(self):
        return None 

###############################
band = BandStructure()