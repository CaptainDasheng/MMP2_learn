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
for i in xrange(0, 1):
    print i
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/NH2CH2COOH.xyz').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/TlCu3Se2.cif').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/TlCu3Se2.vasp').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/In4Te3.cif').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/In4Te3.vasp').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/CONTCAR').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/CONTCAR2', dtype='poscar').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/CONTCAR3', dtype='poscar').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/CONTCAR4', dtype='poscar').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/AgHg3O6Sb.cif').run()
    #raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/AgHg3O6Sb.vasp').run()
    raw=Read('/home/fu/workspace/jump2-0.2.0/examples/structures/K2NaAlF6.vasp').run()
    s1=Structure().create(raw, isPersist=False)
    #s1=Structure().create(raw, isPersist=True, species=['Fe2+', 'Fe3+', 'Fe3+', 'Sb-', 'La2+'])
    print s1
    
    Als=s1.get_element('Al').atoms
    Nas=s1.get_element('Na').atoms
    Ks=s1.get_element('K').atoms
    Fs=s1.get_element('F').atoms
    
    comb=list(itertools.combinations(Als, 3))
    for c in comb:
        s0=StructureFactory(s1).substitute_atoms(c, ['Cs', 'Ba', 'Ir'], isUpdatedInfo=True, isPersist=False).structure
        print s0
    #Write(s1, '/home/fu/workspace/jump2-0.2.0/examples/structures/K2NaAlF6-output.cif', dtype='cif').run()

    s=StructureFactory(s1).del_atoms([1,2,3], isUpdatedInfo=True, isPersist=False).structure
    print s
    
end=get_time()
print '%s: starting...' %start
print '%s: end...' %end

raw2=Read('/home/fu/workspace/Perovskite/StrainedPerovskite/rp/6/0/FAPBI3-rp.vasp').run()
s2=Structure().create(raw2, isPersist=False)
print s2.lattice, s2.get_element('H').atoms
print ''

s=StructureFactory(s2).vacuum([0,0,20,'c'], isUpdatedInfo=True, isPersist=False).structure
Write(s, '/home/fu/workspace/Perovskite/StrainedPerovskite/rp/6/0/FAPBI3-rp-2.vasp', dtype='vasp').run()
print s2.lattice, s2.get_element('H').atoms
print s.lattice, s.get_element('H').atoms
print ''

s_2=StructureFactory(s2).translation(s2.atoms, [0,0,0.1], isUpdatedInfo=True, isPersist=False).structure
Write(s_2, '/home/fu/workspace/Perovskite/StrainedPerovskite/rp/6/0/FAPBI3-rp-3.vasp', dtype='vasp').run()
print s2.lattice, s2.get_element('H').atoms
print s_2.lattice, s_2.get_element('H').atoms

from utils.convert import translation, rotation

position=[11,4,0, 'd']
direction=[0,1,0, 'd']

print 't', translation(position, direction)

axis=[1,1,0, 'd']
theta=[180,'r']
origin=[9,4,0, 'd']
print 'r', rotation(position, axis, theta,origin=origin)
print ''

print s2.get_element('H').atoms
s_3=StructureFactory(s2).rotation([s2.get_element('H').atoms[0]], axis=[1,1,0,'d'], theta=[180,'d'], isUpdatedInfo=True, isPersist=False).structure
print s_3.get_element('H').atoms



raw2=Read('/home/fu/workspace/Perovskite/StrainedPerovskite/rp/6/0/CONTCAR.vasp').run()
s2=Structure().create(raw2, isPersist=False)

del_atoms=[]
for atom in s2.atoms:
    position=atom.position
    if position[-1] < 0.47:
        del_atoms.append(atom)
s=StructureFactory(s2).del_atoms(del_atoms, isUpdatedInfo=True, isPersist=False).structure

add_atoms=[]
for atom in s.atoms:
    position=atom.position
    if position[-1] > 0.52:
        formated_atom=atom.to_formated()
        formated_atom[-2]=1.0-formated_atom[-2]
        add_atoms.append(formated_atom)

s=StructureFactory(s).add_atoms(add_atoms, isUpdatedInfo=True, isPersist=False).structure
#s2.update(isPersist=False)
Write(s, '/home/fu/workspace/Perovskite/StrainedPerovskite/rp/6/0/CONTCAR-2.vasp', dtype='vasp').run()













