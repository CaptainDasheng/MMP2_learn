


# script doc for vasp #

#   vasp_input = """\
#   from jump2set import factory

#   vasp = factory('Vasp')

#   #====parameters==============
#   vasp.prec  = 'Normal'
#   vasp.encut = 300
#   vasp.mesh  = 0.20
#   vasp.vdw   = 'd2'

#   #===sudopotential&vasp=======
#   vasp.pbe = '/home/gordon/workplace/pbe'
#   vasp.program = '/home/gordon/workplce/vasp'

#   #========structure===========
#   from job import Job
#   jobs = Job()

#   for d in  'test':
#       
#       jobs = jobs / 'home/' / d
#       jobs.functional = vasp
#       jobs.tasks = 'scf'
#       jobs.structure =  structure.read('POSCAR')


#   """




def construct_input(functional, filename):


     if functional.lower() == 'vasp':
           with open(filename, 'w') as f:

               f.write(vasp_input)
               f.close()

vasp_input="""\
#!/usr/bin/env python   
from moduleFactory import *
import os 
from utils.read import Read
import cPickle as pickle 
from pool import Pool
#functional to prepare the input for calculation % 
def jump_input():
    
    #three section should be modify at the first calculation:
    #==================================================================
    #section 1: which program you used %

    vasp = Factory.factory('vasp')
    
    #==================================================================
    #section 2: the basic params you used and tasks you would do%

    vasp.energy = 1e-5         # the energy thread %
    vasp.force  = 1e-2         # the force thread if you relax cell %
                               
    vasp.task   = 'scf' 
    #vasp.task   = 'cell volume ions scf band optics' 
                               # any tasks you want to do, 
                               # default is 'scf'.

    #func.vdw = 'b86'
    #func.plusU = 
    #func.magm = 'b86'
    #=================================================================
    #section 3: the initial structure you consider % 
                # you can use the for/while loop to do;
                #here is a simple for loop.
    
    pool = Pool('test.pool')
    pool.functional = vasp

    for d in os.walk('./test').next()[2]: 

        structure =Read('test/'+d, ftype='vasp').getStructure()

        pool.structure =structure

        job = pool / d

    #================================================================
    pool.save_pool()
"""
