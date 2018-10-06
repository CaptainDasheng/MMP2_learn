#!/usr/bin/env python

from jump2.utils.io import runtime
from jump2.functional import Factory 
from jump2.structure.structure import read
from jump2.pool import Pool
import os 

#functional to prepare the input for calculation % 
def jump_input():
    func = Factory.factory('vasp')
    
    #the basic params you used and tasks you would do%
    func.functional = 'pbe'
    func.energy = 1e-4         # the energy thread %
    func.force  = 1e-2         # the force thread if you relax cell %
    func.task   = 'volume ions scf band dos' 
    func.kpoints   = 0.23
    #func.program ='/home/xgzhao/apps/bin/vasp.5.3.5_scf'
    func.program ='/share/apps/vasp5.2/vasp'
    #func.task   = 'cell volume ions scf band optics' 
                               # any tasks you want to do, 
                               # default is 'scf'.
    #func.vdw = 'b86'
    
    #the initial structure you consider % 
    # you can use the for/while loop to do;
    #here is a simple for loop.
    
    pool = Pool()
    print "initial  time:", runtime()
    pool.functional = func
    for d in os.walk('./perovskites_structure').next()[2]:
        structure=read('perovskites_structure/'+d)
        #structure='perovskites_structure/'+d
        #structure =Read('perovskites_structure/'+d, type='poscar').getStructure()
        pool.functional.structure = structure
        job = pool / ('results/'+d.split('_')[-1])

    pool.save('test.dict', overwrite=True)
    print "finished time:", runtime()
