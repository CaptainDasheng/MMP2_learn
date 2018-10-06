#!/usr/bin/env python
#coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import os
import math

class PlotBand():
    def __init__(self):
        self.path = os.getcwd();
        
    # 得到晶胞基矢和倒异基矢
    def getLatticeVector(self):
        infile = open(self.path + "/POSCAR")
        string = infile.readline() # comment line

        string = infile.readline()
        latticeConstant = float(string)
        
        a=[] # lattice vector
        for i in range(0, 3):
            string = infile.readline()
            temp = np.array([float(s0)*latticeConstant for s0 in string.split()])
            if (a == []):
                a = temp
            else:
                a = np.vstack([a, temp])

        volume = a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + \
            a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - \
            a[0][1]*a[1][0]*a[2][2] - a[0][2]*a[1][1]*a[2][0]
        
        b = np.zeros((3,3))
        for i in range(0,3):
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
            c[0] = a[j][1]*a[k][2] - a[j][2]*a[k][1]
            c[1] = a[j][2]*a[k][0] - a[j][0]*a[k][2]
            c[2] = a[j][0]*a[k][1] - a[j][1]*a[k][0]
            for j in range(0, 3):
                b[i][j] = 2*math.pi*c[j]/volume
       
        return a, b 

    # get fermi energy level from OUTCAR
    def getEfermi(self):
        string = os.popen("grep 'E-fermi' OUTCAR").readline()
        Efermi = float(string.split()[2])
        return Efermi

    # get band from EIGENVAL
    # kpoints: n*3 array
    # bands: nbands*nkpoints array
    def getBand(self):

        Efermi = self.getEfermi()
        
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
        kpoints = []
        counter = 0
	maxer=-100.0
	xbar=0.0
        for i in range(0, totalKpoints):
            infile.readline() # blank line
            
            # get kpoint coordinate
            string = infile.readline()
            kp = np.array([float(s0) for s0 in string.split()[:3]])
            if (kpoints == []):
                kpoints = kp
            else:
                kpoints = np.vstack([kpoints, kp])

            # get band value
            bk = [] # bands of a kpoint
            for j in range(0, totalBands):
                string = infile.readline()
                temp = np.array(float(string.split()[1]) - Efermi)
		if temp<=0 and temp>=maxer:
			maxer=temp
                if (bk == []):
                    bk = temp
                else:
                    bk = np.vstack([bk, temp])
                    
            # bk->bands
	    xbar=0.0 - maxer
            if (bands == []):
                bands = bk 
            else:
                bands = np.hstack([bands, bk])
            counter += 1

        if (counter != totalKpoints):
            print "error! read number of kpoints", counter

        return kpoints, bands, xbar

    # claculate the distance of two kpoints
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
	plt.text(float(distances[3]),6,"Two layers", fontsize=20, bbox=dict(facecolor='r', alpha=1, linewidth=0.0) )
        plt.tight_layout()
        plt.savefig("band_two.png", dpi=300)

    # output
    
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
band = PlotBand()
band.plotBand()
