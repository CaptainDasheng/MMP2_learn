'''
Created on Oct 20, 2017

@author: Yuhao Fu
'''
import os
os.environ['DJANGO_SETTINGS_MODULE']='db.settings'
import django
django.setup()

from utils.fetch import get_time

from iostream.read import Read
from iostream.write import Write

from materials.prototype import Prototype
from materials.structure import Structure
from materials.molStructure import MolStructure
from materials.structureFactory import StructureFactory

import pymatgen as mg
import itertools

#si=mg.Element('Si')
#print mg.Composition('FeO3')
#s=mg.Structure.from_file('/home/fu/workspace/jump2-0.2.0/examples/structures/TlCu3Se2.cif')
#s2=mg.Structure.from_file('/home/fu/workspace/jump2-0.2.0/examples/structures/CONTCAR4')
#print s
#print s2
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#finder=SpacegroupAnalyzer(s2)
#print finder.get_hall()


#exit()

start=get_time()
#print '%s: starting...' %start
raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/MoS2.vasp').run()
s1=Structure().create(raw, isPersist=False)

datoms=[]
atoms=s1.atoms
for atom in atoms:
    if atom.position[2] < 0.1 or atom.position[2] > 0.9:
        datoms.append(atom)
print datoms

s=StructureFactory(s1).del_atoms(datoms, isUpdatedInfo=True, isPersist=False).structure
Write(s, '/home/fu/workspace/jump2-0.2.0/examples/structures/MoS2-del.vasp', dtype='vasp').run()
#s=StructureFactory(s1).supercell(dim=[5,5,1], isUpdatedInfo=True, isPersist=False).structure
#s=StructureFactory(s).vacuum([0,0,30,'c'], isUpdatedInfo=True, isPersist=False, isCenter=True).structure
s=StructureFactory(s).vacuum([0,0,30,'c'], isUpdatedInfo=True, isPersist=False, distance=30).structure
Write(s, '/home/fu/workspace/jump2-0.2.0/examples/structures/MoS2-supercell.vasp', dtype='vasp').run()













