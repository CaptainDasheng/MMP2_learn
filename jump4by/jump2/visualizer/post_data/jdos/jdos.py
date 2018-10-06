#!/usr/bin/env python 

import numpy as np


class PlotJDOS(object):

    def __init__(self, ax, path, shift=0.0, color='b'):

        self.path = path
	energy, epsilon = self.get_imaginary()
	energy = np.array(energy) + shift
	self.plot_jdos(ax, energy, epsilon, color)
    
    def get_imaginary(self):
	
	energy = []
        epslon = []
        
        with open(self.path, 'r') as f:
            for i in xrange(0,3): 
                line = f.readline()
            while True:
	        line = f.readline()
                if not line:
                    break
                else:
                    eg = float(line.split()[0])
                    energy.append(eg) 
                    ep = float(line.split()[1])
                    epslon.append(ep)

	return energy, epslon


    def plot_jdos(self, ax, energy, jdos, color):

        ax.plot(energy, jdos, '-', color = color, linewidth=3.0)

