from __future__ import unicode_literals

import math
from fractions import gcd
from itertools import groupby
import numpy as np
import spglib

from django.db import models

from init_elements import InitializeMolElement

from molComposition import MolComposition
from molElement import MolElement
from molAtom import MolAtom


class MolStructureError(Exception):
    pass
    
    
class MolStructure(models.Model,object):
    '''
    molecule structure
    
    Note that:
        coordinate type of positions can only be Cartesian for molecule throughout the interior of the jump2.

    '''
    # relationships
    composition=models.ForeignKey('MolComposition', null=True)
    elements=models.ManyToManyField('MolElement')
    
    label=models.CharField(null=True, blank=True, max_length=80)
    
    natoms=models.IntegerField(blank=True, null=True)
    ntypes=models.IntegerField(blank=True, null=True)
    
    volume=models.FloatField(blank=True, null=True)
    volume_pa=models.FloatField(blank=True, null=True)
    
    
    class Meta:
        app_label='material'
        db_table='molStructure'
        default_related_name='structures'
    
    def __unicode__(self):
        return self.composition.formula
    
    
    def append(self, structure, **kwargs):
        """
        create new molecular structure, append to jump2db
        
        Arguments:
            structure: dictionary of a molecular structure
            
                elements:  ['Ca', 'Fe', 'Sb']
                numbers:   [2, 8, 24]
                positions: [[a1_x,a1_y,a1_z],
                            [a2_x,a2_y,a2_z],
                            [a3_x,a3_y,a3_z],
                                        ...]
                ox (optional):?
                
            kwargs:
                ?
        """
        
        
        # 0. element         
        element_symbols=structure['elements']
        numbers=structure['numbers'] # atomic number of each element
        if len(element_symbols) != len(numbers):
            raise MolStructureError('inconsistent between elements and numbers!')
        
        for e in element_symbols:
            if MolElement.objects.filter(symbol=e).count() == 0: # nonexistent element
                InitializeMolElement(e)
        
        # 2. composition
        # formula
        formula=''
        mass=0
        for i in xrange(0, len(element_symbols)):
            # formula
            if numbers[i] == 1:
                formula += element_symbols[i]
            else:
                formula += element_symbols[i]+str(numbers[i])
            # mass
            mass += MolElement.objects.get(symbol=element_symbols[i]).mass*(numbers[i])
            
        comp, iscreate=MolComposition.objects.get_or_create(formula=formula)
        
        for i in xrange(0, len(element_symbols)):
            comp.elements.add(MolElement.objects.get(symbol=element_symbols[i]))
        comp.mass=mass
        comp.save()
        
        # 3. structure
        #struc, iscreate = Structure.objects.get_or_create(composition=comp)
        struc=MolStructure.objects.create(composition=comp) # Note: can't use get_or_create method, need create new structure object

        for i in xrange(0, len(element_symbols)):
            struc.elements.add(MolElement.objects.get(symbol=element_symbols[i]))
        
        # label       
                
        struc.ntypes=len(element_symbols)
        struc.natoms=np.sum(numbers)

        # volume
        # volume_pa
        
        struc.save()
        
        # 4. atom
        atom_set=[] # atom's object set
        atom_index=0 # index of atoms
        for i in xrange(0, len(element_symbols)): # element
            for j in xrange(0, numbers[i]): # number of element
                #atom, iscreate = Atom.objects.get_or_create(structure=struc,
                atom=MolAtom.objects.create(structure=struc,
                                         element=MolElement.objects.get(symbol=element_symbols[i]),
                                         x=structure['positions'][atom_index][0],
                                         y=structure['positions'][atom_index][1],
                                         z=structure['positions'][atom_index][2]) # Note: coordinate type of positions can only be Cartesian.
                atom_set.append(atom)
                atom_index += 1
                
        return struc
    
    def pop(self):
        """
        remove this molecular structure object from jump2db
        
        Return:
            status of removing (True/False)
        """
        # composition
        comp=self.composition
        comp.structures.remove(self)
        if comp.structures.all().count() == 0:
        #if (comp.structures.count() == 0) and (comp.prototypes.count() == 0): whether need to consider the prototype?
            comp.delete()
        
        # element
        for element in self.elements.all():
            self.elements.remove(element)
            
        # atom
        MolAtom.objects.filter(structure=self).delete()         
                      
        self.delete()
    
    def clone(self):
        """
        duplicate a molecular structure object
        
        Return:
            duplicated molecular structure object
        """
        new=MolStructure()
        new.save()
        
        new.composition=self.composition
        
        # element
        for element in self.elements.all():
            new.elements.add(element)
            
        new.label=self.label
        
        new.natoms=self.natoms
        new.ntypes=self.ntypes
        
        new.volume=self.volume
        new.volume_pa=self.volume_pa
        
        new.save()
        
        # atom
        for atom in self.atoms.all():
            new_atom, iscreate=MolAtom.objects.get_or_create(structure=new,
                                                          element=atom.element,
                                                          ox=atom.ox,
                                                          x=atom.x,
                                                          y=atom.y,
                                                          z=atom.z,
                                                          magmom=atom.magmom,
                                                          charge=atom.charge,
                                                          volume=atom.volume)

        return new
    
    def formatting(self, dtype='xyz', **kwargs):
        """
        translation molecular structure object to a special formatting object. i.e s2xyz s2mol ...
        output a dict array
        
        Arguments:
            dtype: destination type. i.e. xyz or mol
            
            kwargs:
                ?
        """
        if dtype.lower() == 'xyz':
            # elements, numbers and position, constraint(optional)
            elements=[]
            numbers=[]
            positions=[]
            constraints=[]
            for i in xrange(0, self.elements.all().count()):
                elements.append(self.elements.all()[i].symbol)
                numbers.append(self.elements.all()[i].atoms.filter(structure_id=self.id).count())
                for j in xrange(0, numbers[i]):
                    positions.append(self.elements.all()[i].atoms.filter(structure_id=self.id)[j].position)

            elements=np.array(elements)
            numbers=np.array(numbers)
            positions=np.array(positions) 
     
            #poscar=(comment,lattice,elements,numbers,type,positions,constraints)
            molecule={'elements':elements,
                    'numbers':numbers,
                    'type':type,
                    'positions':positions}
            
            return molecule
        
        elif dtype.lower() == 'mol':    
            pass
            
            
    def output(self, path, dtype='xyz'):
        """
        output structure specified structural type to disk
        
        
        Arguments:
            path: path of output. Note that path need to contain the name of file. i.e. /home/xx/POSCAR
            dtype: destination type. i.e. poscar, cif or cell
        """
        output=open(path, 'w')
        
        if dtype.lower() == 'xyz':
            molecule=self.formatting(dtype='xyz')
            
            output.write('%5d\n' %self.natoms)
            output.write(self.composition.formula+'\n')
            
            # element and number
            elements=molecule['elements']
            numbers=molecule['numbers']
            positions=molecule['positions']
            atom_index=0
            for i in xrange(0, elements.shape[0]):
                for j in xrange(0, numbers[i]):
                    output.write('%5s %22.16f %22.16f %22.16f\n' %(elements[i], positions[atom_index][0], positions[atom_index][1], positions[atom_index][2]))
                    atom_index += 1
                    
        elif dtype.lower() == 'mol':
            pass
        
        
        
        
        
        