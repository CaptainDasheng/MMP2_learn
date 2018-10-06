from __future__ import unicode_literals

import math
from fractions import gcd
from itertools import groupby
import numpy as np
import spglib

from django.db import models

from init_elements import InitializeElement

from prototype import Prototype

from composition import Composition
from element import Element
from species import Species
from atom import Atom, cartesian2direct, direct2cartesian
from site import Site

from spacegroup import Translation
from spacegroup import Rotation
from spacegroup import Operation
from spacegroup import Spacegroup
from spacegroup import WyckoffSite

import sys, os
import symbol
pro_dir=os.getcwd()
sys.path.append(pro_dir)
sys.path.append(pro_dir+'/../')
from utils.math import normalizingCoordinate, extractingCoordinate, convertingCoordinate


class StructureError(Exception):
    pass
    
    
class Structure(models.Model,object):
    '''
    crystal structure
    
    Note that:
        type of atomic position is 'Direct' (default) throughout the interior of the jump2.


    '''
    # relationships
    composition=models.ForeignKey('Composition', null=True)
    elements=models.ManyToManyField('Element')
    species=models.ManyToManyField('Species')
    spacegroup=models.ForeignKey('Spacegroup', blank=True, null=True)
    prototype=models.ForeignKey('Prototype', blank=True, null=True, related_name='+')
    
    label=models.CharField(null=True, blank=True, max_length=80)
    comment=models.CharField(default='', max_length=80) # comment line (1st) in POSCAR file. default=composition?
    
    natoms=models.IntegerField(blank=True, null=True)
    nsites=models.IntegerField(blank=True, null=True)
    ntypes=models.IntegerField(blank=True, null=True)
    
    multiple=models.IntegerField(blank=True, null=True) # number of formula
 
    # lattice vector
    x1=models.FloatField(null=True)
    y1=models.FloatField(null=True)
    z1=models.FloatField(null=True)
    x2=models.FloatField(null=True)
    y2=models.FloatField(null=True)
    z2=models.FloatField(null=True)
    x3=models.FloatField(null=True)
    y3=models.FloatField(null=True)
    z3=models.FloatField(null=True)
    
    volume=models.FloatField(null=True)
    volume_pa=models.FloatField(null=True)
    
    energy=models.FloatField(blank=True, null=True) # total energy of structure
    energy_per_formula=models.FloatField(blank=True, null=True)
    energy_per_atom=models.FloatField(blank=True, null=True)
    
    #structure_ptr_id=models.IntegerField(blank=True, null=True) # why need to add the filed. ptr means pointer
    
    class Meta:
        app_label='material'
        db_table='structure'
        default_related_name='structures'
    
    def __unicode__(self):
        return self.composition.formula
    
    
    def append(self, structure, **kwargs):
        """
        create new structure, append to jump2db
        Note: the value of atomic coordinate (Direct) is between 0.0 and 1.0 by default.
        
        Arguments:
            structure: dictionary of a structure
            
            lattice: [[x1,y1,z1],
                     [x2,y2,z2],
                     [x3,y3,z3]]
            elements:  ['Ca', 'Fe', 'Sb']
            numbers:   [2, 8, 24]
            positions: [[a1_x,a1_y,a1_z],
                        [a2_x,a2_y,a2_z],
                        [a3_x,a3_y,a3_z],
                                    ...]
            ox (optional):?
        
            kwargs:
                isNormalizingCoordinate (default=True): consider the translation periodicity when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                # symmetrized parameter
                symprec: default 1e-5
                angle_tolerance: a experimental argument that controls angle tolerance between basis vectors. the default value is -1.0.
                #hall_number: default 0
                
                energy: total energy of structure
                energy_per_formula: energy per formula
                energy_per_atom: energy per atom
           
        Return:
            created structure object
           
        """
        # step:
        # 0. element
        #        z, name, symbol, group, period, mass, density(?), volume(?)
        # 1. species (optional)(?)
        #        ->element, name, ox
        # 2. composition
        #        ->element, formlar, generic(?), mass
        # 3. structure
        #        ->composition, ->element, ->species(?), spacegroup(step 5), ->prototype(?), label, comment, natoms, nsites(step 6), ntypes, x1..z3, volume, volume_pa
        # 4. atom
        #        ->structure, ->element, ->site(step 6), ox(?), x..z, fx..fz(?), magmom(?), occupancy(?), wyckoff(step 6), cx..cz, species(undecided parameter)(?)
        # 5. spacegroup; translation; rotation; operation
        #    cell, calling spglib
        #    spacegroup
        #        ->operations(step 5), number, hm(?), hall, pearson(?), schoenflies(?), lattice_system(?), centerosymmetric(?)
        #    translation
        #        x..z
        #    rotation
        #        a11..a33
        #    operation
        #        ->rotation, ->translation
        #    spacegroup.operations, structures.spacegroup
        # 6. wyckoffSite; site
        #    wyckoffSite
        #        ->spacegroup, symbol, multiplicity, x..z
        #    site
        #        ->structure, ->wyckoff, x..z, coordiantion_number(?)
        #    structure.nsites, atom.wyckoff, atom.site 

        # 0. element         
        element_symbols=structure['elements']
        numbers=structure['numbers'] # atomic number of each element
        if len(element_symbols) != len(numbers):
            raise StructureError('inconsistent between elements and numbers!')
        
        multi=reduce(gcd, numbers) # number of formula
        
        for e in element_symbols:
            if Element.objects.filter(symbol=e).count() == 0: # nonexistent element
                InitializeElement(e)

        # 1. species(optional)
        
        # 2. composition
        # formula
        formula=''
        mass=0
        for i in xrange(0, len(element_symbols)):
            # formula
            if numbers[i]/multi == 1:
                formula += element_symbols[i]
            else:
                formula += element_symbols[i]+str(numbers[i]/multi)
            # mass
            mass += Element.objects.get(symbol=element_symbols[i]).mass*(numbers[i]/multi)
            
        comp, iscreate=Composition.objects.get_or_create(formula=formula)
        
        for i in xrange(0, len(element_symbols)):
            comp.elements.add(Element.objects.get(symbol=element_symbols[i]))
        comp.mass=mass
        comp.save()
        
        # 3. structure
        #struc, iscreate = Structure.objects.get_or_create(composition=comp)
        struc=Structure.objects.create(composition=comp) # Note: can't use get_or_create method, need create new structure object

        for i in xrange(0, len(element_symbols)):
            struc.elements.add(Element.objects.get(symbol=element_symbols[i]))
            
        # species
        
        # prototype
        
        # label
        
        if 'comment' in structure:
            struc.comment=structure['comment']        
                
        struc.ntypes=len(element_symbols)
        struc.natoms=np.sum(numbers)

        struc.x1=structure['lattice'][0][0]
        struc.y1=structure['lattice'][0][1]
        struc.z1=structure['lattice'][0][2]
        struc.x2=structure['lattice'][1][0]
        struc.y2=structure['lattice'][1][1]
        struc.z2=structure['lattice'][1][2]
        struc.x3=structure['lattice'][2][0]
        struc.y3=structure['lattice'][2][1]
        struc.z3=structure['lattice'][2][2]
        struc.volume=struc.calculatedVolume()
        struc.volume_pa=struc.volume/struc.natoms
        
        struc.multiple=multi
        
        if 'energy' in kwargs:
            struc.energy=float(kwargs['energy'])
            struc.energy_per_formula=struc.energy/struc.multiple
            struc.energy_per_atom=struc.energy/struc.natoms
            
        if 'energy_per_formula' in kwargs:
            struc.energy_per_formula=float(kwargs['energy_per_formula'])
            struc.energy=struc.energy_per_formula*struc.multiple
            struc.energy_per_atom=struc.energy/struc.natoms
            
        if 'energy_per_atom' in kwargs:
            struc.energy_per_atom=float(kwargs['energy_per_atom'])
            struc.energy=struc.energy_per_atom*struc.natoms
            struc.energy_per_formula=struc.energy/struc.multiple
            
        struc.save()
        
        # 4. atom
        # translation periodicity
        isNormalizingCoordinate=True
        if 'isNormalizingCoordinate' in kwargs:
            isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
        
        atom_set=[] # atom's object set
        atom_index=0 # index of atoms
        if structure['type'] == 'Direct':
            for i in xrange(0, len(element_symbols)): # element
                for j in xrange(0, numbers[i]): # number of element
                    position=structure['positions'][atom_index]
                    if isNormalizingCoordinate: # remove translation periodicity
                        position=normalizingCoordinate(position)
                        atom=Atom.objects.create(structure=struc,
                                                 element=Element.objects.get(symbol=element_symbols[i]),
                                                 x=position[0],
                                                 y=position[1],
                                                 z=position[2])
                    else:
                        atom=Atom.objects.create(structure=struc,
                                                 element=Element.objects.get(symbol=element_symbols[i]),
                                                 x=position[0],
                                                 y=position[1],
                                                 z=position[2])
                    atom_set.append(atom)
                    atom_index += 1
        elif structure['type'] == 'Cartesian':
            for i in xrange(0, len(element_symbols)): # element
                for j in xrange(0, numbers[i]): # number of element
                    position=cartesian2direct(struc, structure['positions'][atom_index])
                    if isNormalizingCoordinate: # remove translation periodicity
                        position=normalizingCoordinate(position)
                        atom=Atom.objects.create(structure=struc,
                                                 element=Element.objects.get(symbol=element_symbols[i]),
                                                 x=position[0],
                                                 y=position[1],
                                                 z=position[2])
                    else:
                        atom=Atom.objects.create(structure=struc,
                                                 element=Element.objects.get(symbol=element_symbols[i]),
                                                 x=position[0],
                                                 y=position[1],
                                                 z=position[2])
                    atom_set.append(atom)            
                    atom_index += 1
        else:
            raise StructureError('unknown value in type(Direct/Cartesian)!')
                            
        # ox
        
        # fx, fy, fz
        
        # magmom
        
        # occupancy
        
        # cx, cy, cz
        if ('constraints' in structure.keys()) and (structure['constraints'] != []):
            for i in xrange(0, len(atom_set)): # atom
                atom_set[i].cx=structure['constraints'][i][0]
                atom_set[i].cy=structure['constraints'][i][1]
                atom_set[i].cz=structure['constraints'][i][2]
                atom_set[i].save()
        
        # species
        
        # 5. spacegroup; translation; rotation; operation
        # cell for input parameter of spglib
        cell=struc.formatting('cell')
        cell=(cell['lattice'], cell['positions'], cell['numbers'])
        
        # default
        stmp=1e-5 # symprec
        atmp=-1.0 # angle_tolerance
        #hall_number = 0
        if 'symprec' in kwargs:
            stmp=kwargs['symprec']
        if 'angle_tolerance' in kwargs:
            atmp=kwargs['angle_tolerance']
        #if 'hall_number' in kwargs:
        #    hall_number=kwargs['hall_number']
            
        dataset=spglib.get_symmetry_dataset(cell, symprec=stmp, angle_tolerance=atmp)

        spacegroup, iscreate=Spacegroup.objects.get_or_create(number=dataset["number"],
                                                              hall=dataset['hall'])
        # hm
        # pearson
        # schoenflies
        # lattice_system
        # centerosymmetric
        
        for i in xrange(0, len(dataset['translations'])):
            translation, iscreate=Translation.objects.get_or_create(x=dataset["translations"][i][0],
                                                                    y=dataset["translations"][i][1],
                                                                    z=dataset["translations"][i][2])
            rotation, iscrate=Rotation.objects.get_or_create(a11=dataset["rotations"][i][0][0],
                                                             a12=dataset["rotations"][i][0][1],
                                                             a13=dataset["rotations"][i][0][2],
                                                             a21=dataset["rotations"][i][1][0],
                                                             a22=dataset["rotations"][i][1][1],
                                                             a23=dataset["rotations"][i][1][2],
                                                             a31=dataset["rotations"][i][2][0],
                                                             a32=dataset["rotations"][i][2][1],
                                                             a33=dataset["rotations"][i][2][2])
            operation, iscreate=Operation.objects.get_or_create(rotation=rotation,
                                                                translation=translation)
            spacegroup.operations.add(operation)
        
        # structure.spacegroup
        struc.spacegroup=spacegroup
        struc.save()

        # 6. wyckoffSite; site
        wyckoff=dataset["wyckoffs"]
        equivalent_atoms=dataset['equivalent_atoms']
        wyckoff_symbol=[]
        wyckoff_symbol_number=[]
        
        for key, group in groupby(equivalent_atoms):
            wyckoff_symbol.append(wyckoff[key])
            wyckoff_symbol_number.append(len(list(group)))

        site_index=[]
        for i in xrange(0, len(wyckoff_symbol_number)):
            site_index.append(sum(wyckoff_symbol_number[:i]))

        wyckoff_set=[] # wyckoff's object set
        site_set=[] # site's object set
        for i in xrange(0, len(wyckoff_symbol)):
            wyck, iscreate=WyckoffSite.objects.get_or_create(spacegroup=spacegroup,
                                                             symbol=wyckoff_symbol[i],
                                                             multiplicity=wyckoff_symbol_number[i],
                                                             x=structure['positions'][site_index[i]][0],
                                                             y=structure['positions'][site_index[i]][1],
                                                             z=structure['positions'][site_index[i]][2])
            wyckoff_set.append(wyck)
            #site, iscreate = Site.objects.get_or_create(structure=struc, 
            site=Site.objects.create(structure=struc,
                                     wyckoff=wyck,
                                     x=structure['positions'][site_index[i]][0],
                                     y=structure['positions'][site_index[i]][1],
                                     z=structure['positions'][site_index[i]][2])
            
            # coordiantion_number
            
            site_set.append(site)
        
            
        # structure.nsites, atom.wyckoff, atom.site 
        struc.nsites=len(site_set)
        struc.save()
          
        atom_index=0 # index of atoms
        for i in xrange(0, len(wyckoff_symbol_number)):
            for j in xrange(0, wyckoff_symbol_number[i]):
                atom_set[atom_index].wyckoff=wyckoff_set[i]
                atom_set[atom_index].site=site_set[i]
                atom_set[atom_index].save()
                atom_index += 1
                
        return struc

    def pop(self, **kwargs):
        """
        remove this structure object from jump2db
        
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
            
        # species
        for species in self.species.all():
            self.species.remove(species)
        
        # site
        Site.objects.filter(structure=self).delete()
        
        # atom
        Atom.objects.filter(structure=self).delete()         
            
        # spacegroup
        spacegroup=self.spacegroup
        spacegroup.structures.remove(self)
        if spacegroup.structures.all().count() == 0:
            # operation
            for operation in spacegroup.operations.all():
                if operation.spacegroups.all().count() == 1:
                    if (operation.rotation.operations.all().count() == 1) and (operation.translation.operations.all().count() == 1):
                        operation.rotation.delete()
                        operation.translation.delete()
                    elif operation.rotation.operations.all().count() == 1:
                        operation.rotation.delete()
                    elif operation.translation.operations.all().count() == 1:
                        operation.translation.delete()
                    operation.delete()
                
            # wyckoffSite
            for wyckoff in spacegroup.wyckoffSites.all():
                if wyckoff.atoms.all().count() == 0:
                    wyckoff.delete()
                
            spacegroup.delete()
            
        # prototype
        prototype=self.prototype
        if prototype != None:
            prototype.structures.remove(self)
                        
        self.delete()
    
    def clone(self):
        """
        duplicate a structure object
        
        Return:
            duplicated structure object
        """
        new=Structure()
        new.save()
        
        new.composition=self.composition
        
        # element
        for element in self.elements.all():
            new.elements.add(element)
            
        # species
        for species in self.species.all():
            new.species.add(species)
            
        new.spacegroup=self.spacegroup
        new.prototype=self.prototype
        
        new.label=self.label
        new.comment=self.comment
        
        new.natoms=self.natoms
        new.nsites=self.nsites
        new.ntypes=self.ntypes
        
        # method 1
        #new.x1=self.x1
        #new.x2=self.x2
        #new.x3=self.x3
        #new.y1=self.y1
        #new.y2=self.y2
        #new.y3=self.y3
        #new.z1=self.z1
        #new.z2=self.z2
        #new.z3=self.z3
        # method 2
        new.lattice=self.lattice

        new.volume=self.volume
        new.volume_pa=self.volume_pa
        
        new.save()
        
        # site
        map={} # {old:new,...} old site map to new site.
        for site in self.sites.all():
            new_site, iscreate=Site.objects.get_or_create(structure=new,
                                                          wyckoff=site.wyckoff,
                                                          x=site.x,
                                                          y=site.y,
                                                          z=site.z,
                                                          coordiantion_number=site.coordiantion_number)
            map[site.id]=new_site

        # atom
        for atom in self.atoms.all():
            new_site=map[atom.site.id]
            new_atom, iscreate=Atom.objects.get_or_create(structure=new,
                                                          site=new_site,
                                                          element=atom.element,
                                                          species=atom.species,
                                                          ox=atom.ox,
                                                          x=atom.x,
                                                          y=atom.y,
                                                          z=atom.z,
                                                          fx=atom.fx,
                                                          fy=atom.fy,
                                                          fz=atom.fz,
                                                          cx=atom.cx,
                                                          cy=atom.cy,
                                                          cz=atom.cz,
                                                          magmom=atom.magmom,
                                                          charge=atom.charge,
                                                          volume=atom.volume,
                                                          occupancy=atom.occupancy,
                                                          wyckoff=atom.wyckoff,
                                                          )

        return new

    def minimizeClone(self):
        """
        duplicate minimum information of structure, only including lattice and atoms.
        
        Return:
            duplicated structure object
        """
        new=Structure()
        new.save()
        
        new.composition=self.composition
        
        # element
        for element in self.elements.all():
            new.elements.add(element)
            
        # species
        #for species in self.species.all():
        #    new.species.add(species)
            
        #new.spacegroup=self.spacegroup
        #new.prototype=self.prototype
        
        #new.label=self.label
        #new.comment=self.comment
        
        new.natoms=self.natoms
        #new.nsites=self.nsites
        new.ntypes=self.ntypes
        
        # method 1
        #new.x1=self.x1
        #new.x2=self.x2
        #new.x3=self.x3
        #new.y1=self.y1
        #new.y2=self.y2
        #new.y3=self.y3
        #new.z1=self.z1
        #new.z2=self.z2
        #new.z3=self.z3
        # method 2
        new.lattice=self.lattice

        new.volume=self.volume
        new.volume_pa=self.volume_pa
        
        new.save()
        
        # site
        #map={} # {old:new,...} old site map to new site.
        #for site in self.sites.all():
        #    new_site, iscreate=Site.objects.get_or_create(structure=new,
        #                                                  wyckoff=site.wyckoff,
        #                                                  x=site.x,
        #                                                  y=site.y,
        #                                                  z=site.z,
        #                                                  coordiantion_number=site.coordiantion_number)
        #    map[site.id]=new_site

        # atom
        for atom in self.atoms.all():
            new_site=map[atom.site.id]
            new_atom, iscreate=Atom.objects.get_or_create(structure=new,
                                                          element=atom.element,
                                                          site=new_site,
                                                          ox=atom.ox,
                                                          x=atom.x,
                                                          y=atom.y,
                                                          z=atom.z,
                                                          fx=atom.fx,
                                                          fy=atom.fy,
                                                          fz=atom.fz,
                                                          magmom=atom.magmom,
                                                          occupancy=atom.occupancy,
                                                          wyckoff=atom.wyckoff,
                                                          cx=atom.cx,
                                                          cy=atom.cy,
                                                          cz=atom.cz)

        return new
        
    def formatting(self, dtype='poscar', **kwargs):
        """
        translation structure object to a special formatting object. i.e s2cif s2poscar s2cell...
        output a dict array
        
        Arguments:
            dtype: destination type. i.e. poscar, cif or cell
            
            **kwargs:
                coord_type:type of atomic coordinate (Direc/Cartesian).
                isConstraint:True/False
                
        Return:
            a array of this tidy structure (type:poscar/cif/cell)
            
            poscar: {'comment':comment,
                     'lattice':lattice,
                     'elements':elements,
                     'numbers':numbers,
                     'type':type,
                     'positions':positions,
                     'constraints':constraints(optional)}
                    
            cif(?)
            
            cell: {'lattice':lattice,
                   'positions':positions,
                   'numbers':numbers,
                   'magmoms':magmoms(optional)}
        
        poscar:
            comment: comment of the first line
            lattice=[[x1,y1,z1],
                     [x2,y2,z2],
                     [x3,y3,z3]]
            elements=['Ca', 'Fe', 'Sb']
            numbers=[2, 8, 24]
            type= Direct or Cartesian
            positions=[[a1_x,a1_y,a1_z],
                      [a2_x,a2_y,a2_z],
                      [a3_x,a3_y,a3_z],
                      ...]
            constraints=[[T,T,T], # Selective dynamics (optional)
                        [F,F,F],
                        [T,F,T],
                        ...]
            
        cell:
            cell=(lattice, positions, numbers) or (lattice, positions, numbers, magmoms)
        
            lattice=[[x1,y1,z1],
                     [x2,y2,z2],
                     [x3,y3,z3]]
            positions=[[a1_1,a1_2,a1_3],
                       [a2_1,a2_2,a2_3],
                       [a3_1,a3_2,a3_3],
                       ...]
            numbers=[n_1,n_2,n_3,...] numbers to distinguish atomic species by distinguishing different elements?
            magmoms=[m_1,m_2,m_3,...] (optional)
        """
        if dtype.lower() == 'poscar':
            
            poscar=()
             
            comment=self.comment
            lattice=self.lattice
            
            type='Direct'
            if 'coord_type' in kwargs:
                if kwargs['coord_type'].lower().startswith('c'):
                    type='Cartesian'
                elif kwargs['coord_type'].lower().startswith('d'):
                    type='Direct'
                else:
                    raise StructureError('unknown value in type(Direct/Cartesian)!')
            
            # elements, numbers and position, constraint(optional)
            elements=[]
            numbers=[]
            positions=[]
            constraints=[]
            for i in xrange(0, self.elements.all().count()):
                elements.append(self.elements.all()[i].symbol)
                numbers.append(self.elements.all()[i].atoms.filter(structure_id=self.id).count())
                for j in xrange(0, numbers[i]):
                    if type == 'Direct':
                        positions.append(self.elements.all()[i].atoms.filter(structure_id=self.id)[j].position)
                    else:
                        positions.append(self.elements.all()[i].atoms.filter(structure_id=self.id)[j].direct2cartesian())
                    
                    if 'isConstraint' in kwargs and kwargs['isConstraint']:
                        constraints.append(self.elements.all()[i].atoms.filter(structure_id=self.id)[j].constraint)
            elements=np.array(elements)
            numbers=np.array(numbers)
            positions=np.array(positions) 
            constraints=np.array(constraints)
     
            #poscar=(comment,lattice,elements,numbers,type,positions,constraints)
            poscar={'comment':comment,
                    'lattice':lattice,
                    'elements':elements,
                    'numbers':numbers,
                    'type':type,
                    'positions':positions,
                    'constraints':constraints}
            
            return poscar
             
        elif dtype.lower() == 'cif':
            pass
        elif dtype.lower() == 'cell':
            numbers=[]
            positions=[]
            for i in xrange(0, self.elements.all().count()):
                for j in xrange(0, self.elements.all()[i].atoms.filter(structure_id=self.id).count()):
                    numbers.append(self.elements.all()[i].z)
                    positions.append(self.elements.all()[i].atoms.filter(structure_id=self.id)[j].position)
                    
            numbers=np.array(numbers)
            positions=np.array(positions)
            
            cell=()
            #cell=(self.lattice,
            #      np.array([x.position for x in self.atoms.all()]),
            #      numbers)
            cell={'lattice':self.lattice,
                  'positions':positions,
                  'numbers':numbers}            
            return cell
        else:
            raise StructureError('unknown value in dtype(poscar/cif/cell)!')

    def output(self, path, dtype='poscar', **kwargs):
        """
        output structure specified structural type to disk
        
        
        Arguments:
            path: path of output. Note that path need to contain the name of file. i.e. /home/xx/POSCAR
            dtype: destination type. i.e. poscar, cif or cell
            
            **kwargs:
                coord_type: Direct/Cartesian
                selective_dynamics: boolean value (True/False).
                ##isConstraint: constraint of atom (True/False). <- deprecated
        """
        output=open(path, 'w')
        
        if dtype.lower() == 'poscar':
            poscar={}
            if ('coord_type' in kwargs) and (('selective_dynamics' in kwargs) and kwargs['selective_dynamics']):
                poscar=self.formatting(dtype='poscar', coord_type=kwargs['coord_type'], isConstraint=True)
            elif 'coord_type' in kwargs:
                poscar=self.formatting(dtype='poscar', coord_type=kwargs['coord_type'])
            elif ('selective_dynamics' in kwargs) and kwargs['selective_dynamics']:
                poscar=self.formatting(isConstraint=True)
            else:
                poscar=self.formatting(dtype='poscar')
            
            output.write(poscar['comment']+'\n')
            
            output.write('   1.0'+'\n') # scale value. Note that: don't modify it. If you do it, the lattice will be change.
            
            # lattice
            lattice=poscar['lattice']
            for i in xrange(0, lattice.shape[0]):
                output.write(' %22.16f %22.16f %22.16f\n' %(lattice[i][0], lattice[i][1], lattice[i][2]))
                
            # element
            elements=poscar['elements']
            for i in xrange(0, elements.shape[0]):
                output.write('%5s' %elements[i])
            output.write('\n')
            # number
            numbers=poscar['numbers']
            for i in xrange(0, numbers.shape[0]):
                output.write('%5d' %numbers[i])
            output.write('\n')
                
            # Selective dynamics
            if ('selective_dynamics' in kwargs) and kwargs['selective_dynamics']:
                output.write('Selective dynamics\n')
                
                # type
                output.write(poscar['type']+'\n')
            
                # postions
                positions=poscar['positions']
                constraints=poscar['constraints']
                for i in xrange(0, positions.shape[0]):
                    output.write('%22.16f %22.16f %22.16f' %(positions[i][0], positions[i][1], positions[i][2]))
                    constraint=['T' if s0 else 'F' for s0 in constraints[i]]
                    output.write('%4s %4s %4s\n' %(constraint[0], constraint[1], constraint[2]))
            else:
                # type
                output.write(poscar['type']+'\n')
            
                # postions
                positions=poscar['positions']
                #constraints=poscar['constraints']
                for i in xrange(0, positions.shape[0]):
                    output.write('%22.16f %22.16f %22.16f\n' %(positions[i][0], positions[i][1], positions[i][2]))
                    
        elif dtype.lower() == 'cif':
            pass
        
        elif dtype.lower() == 'cell':
            pass
        else:
            raise StructureError('unknown value in dtype(poscar/cif/cell)!')                
            
        output.close()


    def setEnergy(self, e):
        """
        set the value of energy(total energy of structure), and synchronize the value of energy_per_formula and energy_per_atom
        
        Arguments:
            e: total energy of structure
        """
        e=float(e)
        self.energy=e
        self.energy_per_formula=e/self.multiple
        self.energy_per_atom=e/self.natoms
        
        self.save()
  
    def setEnergyPerFormula(self, e):
        """
        set the value of energy_per_formula, and synchronize the values of energy and energy_per_atom
        
        Arguments:
            e: energy per formula
        """
        e=float(e)
        self.energy_per_formula=e
        self.energy=e*self.multiple
        self.energy_per_atom=self.energy/self.natoms
        
        self.save()
    
    def setEnergyPerAtom(self, e):
        """
        set the value of energy_per_atom, and synchronize the value of energy and energy_per_formula
        
        Arguments:
            e: energy per atom
        """
        e=float(e)
        self.energy_per_atom=e
        self.energy=e*self.natoms
        self.energy_per_formula=self.energy/self.multiple
        
        self.save()

        
    @property
    def lattice(self):
        """
        lattice vector
        
        Return:
            a lattice 2D-array of this stucture. [[x1,y1,z1],
                                                  [x2,y2,z2],
                                                  [x3,y3,z3]]
            
        """
        self._lattice=np.array([[self.x1, self.y1, self.z1],
                         [self.x2, self.y2, self.z2],
                         [self.x3, self.y3, self.z3]])
        
        return self._lattice

    @lattice.setter
    def lattice(self, lattice):
        """
        Arguments:
            lattice:a 3x3 array lattice vector. i.e. [[x1,y1,z1],
                                                      [x2,y2,z2],
                                                      [x3,y3,z3]]
        """
        # check lattice
        lattice=np.array(lattice)
        if lattice.shape != (3,3):
            raise StructureError('invalid lattice!')
        
        self.x1=lattice[0][0]
        self.y1=lattice[0][1]
        self.z1=lattice[0][2]
        self.x2=lattice[1][0]
        self.y2=lattice[1][1]
        self.z2=lattice[1][2]
        self.x3=lattice[2][0]
        self.y3=lattice[2][1]
        self.z3=lattice[2][2]
        self.save()
    
    @property
    def lattice_parameters(self):
        """
        gain the lattice parameters.
        
        Return:
            a array of lattice parameters. [a, b, c, alpha, beta, gamma]
        """
        a=np.linalg.norm(self.lattice[0])
        b=np.linalg.norm(self.lattice[1])
        c=np.linalg.norm(self.lattice[2])
        alpha=np.degrees(np.arccos(np.clip(np.dot(self.lattice[1]/b, self.lattice[2]/c), -1, 1)))
        beta=np.degrees(np.arccos(np.clip(np.dot(self.lattice[0]/a, self.lattice[2]/c), -1, 1)))
        gamma=np.degrees(np.arccos(np.clip(np.dot(self.lattice[0]/a, self.lattice[1]/b), -1, 1)))
        
        return np.array([a, b, c, alpha, beta, gamma])
        
    @property    
    def reciprocal_lattice(self):
        """
        reciprocal lattice vector
        
        Note that: test if the calculated reciprocal lattice is right. 2.math.pi?
        
        Return:
            a reciprocal lattice 2D-array of this structure. [[x1,y1,z1],
                                                              [x2,y2,z2],
                                                              [x3,y3,z3]]
        """
        self._reciprocal_lattice=np.zeros((3,3))
        for i in xrange(0,3):
            if i == 0:
                j=1
                k=2
            elif i ==1:
                j=2
                k=0
            else:
                j=0
                k=1
            tmp=np.zeros(3)
            tmp[0]=self.lattice[j][1]*self.lattice[k][2]-self.lattice[j][2]*self.lattice[k][1]
            tmp[1]=self.lattice[j][2]*self.lattice[k][0]-self.lattice[j][0]*self.lattice[k][2]
            tmp[2]=self.lattice[j][0]*self.lattice[k][1]-self.lattice[j][1]*self.lattice[k][0]
            for j in xrange(0,3):
                self._reciprocal_lattice[i][j]=2*math.pi*tmp[j]/self.volume
                
        return self._reciprocal_lattice
               
    def calculatedVolume(self):
        """
        get the volume of this structure by calculating its lattice vector.
        
        Return:
            calculated volume (float).
        """
        
        lattice=self.lattice
        
        # method 1
        #volume=lattice[0][0]*lattice[1][1]*lattice[2][2] + \
        #    lattice[0][1]*lattice[1][2]*lattice[2][0] + \
        #    lattice[0][2]*lattice[1][0]*lattice[2][1] - \
        #    lattice[0][0]*lattice[1][2]*lattice[2][1] - \
        #    lattice[0][1]*lattice[1][0]*lattice[2][2] - \
        #    lattice[0][2]*lattice[1][1]*lattice[2][0]
            
        # method 2
        volume=np.linalg.det(lattice)
        return volume

    
    def symmetrized(self,symprec=1e-5, angle_tolerance=-1.0):
        """
        find symmetry
        
        Arguments:
            symprec:
            angle_tolerance: a experimental argument that controls angle tolerance between basis vectors.
        
        Return:
            a string of symmetry information. i.e. 'Im3 (204)'
        """
        cell=self.formatting('cell')
        cell=(cell['lattice'], cell['positions'], cell['numbers'])
        dataset=spglib.get_symmetry_dataset(cell)
        
        return '%s (%i)' %(dataset['international'], dataset['number'])
    
    def compare(self, other):
        """
        compare with other structure
        
        
        """
        pass
    
    # override python's method
    def __contains__(self, key, tolerance=1e-3, **kwargs):
        """
        override python's 'in' operator.
        
        Arguments:
            key:
            i.e.
            for element
                struc1.exits('Na') -> element object or None 
            for atom
                [element_symbol, x, y, z, coord_type (optional)] in struc1 -> True/False
            for position
                [x, y, z, coord_type (optional)] in struc1 -> True/False
                
            kwargs:
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                isConstraint (default=False):True/False.
        """
        if isinstance(key, basestring): # element
            return self.__isExistElement(key)
        elif isinstance(key, list):
            if not(isinstance(key[0], int) or isinstance(key[0], float)): # atom
                isNormalizingCoordinate=True
                if 'isNormalizingCoordinate' in kwargs:
                    isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
                isConstraint=False
                if 'isConstraint' in kwargs:
                    isConstraint=kwargs['isConstraint']
                    
                return self.__isExistAtom(key, tolerance, isNormalizingCoordinate=isNormalizingCoordinate, isConstraint=isConstraint)
            else: # site
                isNormalizingCoordinate=True
                if 'isNormalizingCoordinate' in kwargs:
                    isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
                isConstraint=False
                
                return self.__isExistSite(key, tolerance, isNormalizingCoordinate=isNormalizingCoordinate)
        else:
            raise StructureError('unknown key')
        
    def exists(self, key, tolerance=1e-3, **kwargs):
        """
        determine whether the atom or site is in the structure
        
        Arguments:
            key:
            i.e.
            for element
                struc1.exits('Na') -> element object or None
            for atom
                struc1.exits([element_symbol, x, y, z, coord_type (optional)]) -> atomic object or None
            for position
                struc1.exits([x, y, z, coord_type (optional)]) -> atomic object on this site or None
                
            kwargs:
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                isConstraint (default=False):True/False.
        """
        if isinstance(key, basestring): # element
            return self.__isExistElement(key)
        elif isinstance(key, list):
            if not(isinstance(key[0], int) or isinstance(key[0], float)): # atom
                isNormalizingCoordinate=True
                if 'isNormalizingCoordinate' in kwargs:
                    isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
                isConstraint=False
                if 'isConstraint' in kwargs:
                    isConstraint=kwargs['isConstraint']
                    
                return self.__isExistAtom(key, tolerance, isNormalizingCoordinate=isNormalizingCoordinate, isConstraint=isConstraint)
            else: # site
                isNormalizingCoordinate=True
                if 'isNormalizingCoordinate' in kwargs:
                    isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
                isConstraint=False
                
                return self.__isExistSite(key, tolerance, isNormalizingCoordinate=isNormalizingCoordinate)
        else:
            raise StructureError('unknown key')
        
    def __isExistElement(self, element):
        """
        determines whether the element is in the structure
        
        arguments:
            element: symbol of element. i.e. 'Na'
            
        Return:
            return element object.
        """
        if self.elements.filter(symbol=element).count() == 1:
            return self.elements.filter(symbol=element)[0]
        else:
            return None
            
    
    def __isExistAtom(self, atom, tolerance=1e-3, **kwargs):
        """
        determines whether the atom is in the structure.
        
        Arguments:
            atom: a array contained atomic information [element_symbol, x, y, z, coord_type (optional)]. i.e. ['Na', 0.5, 0.5, 0.5](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian'](Cartesian)
            tolerance (default=1e-3): distance tolerance between the position of the input and position of structure, unite is angstrom.
            
            kwargs:
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                isConstraint (default=False):True/False.
                
        Return:
            Return atom object.
        """
        # check data
        atom_set=[]
        element_symbol=atom[0]
        # method 1
        #if not(self.elements.filter(symbol=element_symbol).exists()):
        #    return None
        #else:
        #    atom_set=self.elements.filter(symbol=element_symbol)[0].atoms.filter(structure=self)
        #    if atom_set.count() == 0:
        #        raise StructureError('atomic number of this element equal zero!')
        # method 2
        if self.__isExistElement(element_symbol) != None:
            atom_set=self.elements.filter(symbol=element_symbol)[0].atoms.filter(structure=self)
            
        position=extractingCoordinate(self, atom)
        
        # constraint
        constrained=False # whether constrain atom
        if 'isConstraint' in kwargs:
            constrained=kwargs['isConstraint']
            
        counter=0 # atomic counter satisfying the filter condition
        index=0 # atomic index which need to delete
        tmpa=None # temporary atoms
        for i in xrange(0, atom_set.count()):
            tmpd=position-atom_set[i].position # distance
            
            # remove translation periodicity
            isNormalizingCoordinate=True
            if 'isNormalizingCoordinate' in kwargs:
                isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
            if isNormalizingCoordinate:
                tmpd=normalizingCoordinate(tmpd, range='center')
                    
            distance=np.linalg.norm(tmpd)

            if distance <= tolerance:
                index=i
                tmpa=atom_set[i]
                counter += 1
        try:        
            if counter == 1:
                return tmpa
            elif counter == 0:
                return None
            else:
                raise StructureError('more than one atom satisfies the deletion condition')
        except StructureError:
            raise StructureError()

    def __isExistSite(self, position, tolerance=1e-3, **kwargs):
        """
        determines whether the atom is on this site. Note that the site is all lattice sites, not only the irreducible sites.
        
        Arguments:
            position: a array of atomic position [x, y, z, coord_type (optional)]. i.e. [0.5, 0.5, 0.5](Direc) or [1.3, 3.5, 2.0, 'Cartesian'](Cartesian)
            tolerance (default=1e-3): distance tolerance between the position of the input and position of structure, unite is angstrom.
            
            kwargs:
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                
        Return:
            Return atom object.
        """
        position=convertingCoordinate(self, position)
        
        atom_set=self.atoms.all()
        counter=0 # atomic counter satisfying the filter condition
        index=0 # atomic index which need to delete
        tmpp=None # temporary position
        for i in xrange(0, atom_set.count()):
            tmpd=position-atom_set[i].position # distance
            
            # remove translation periodicity
            isNormalizingCoordinate=True
            if 'isNormalizingCoordinate' in kwargs:
                isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
            if isNormalizingCoordinate:
                tmpd=normalizingCoordinate(tmpd, range='center')
                    
            distance=np.linalg.norm(tmpd)

            if distance <= tolerance:
                index=i
                tmpp=atom_set[i]
                counter += 1
        try:        
            if counter == 1:
                return tmpp
            elif counter == 0:
                return None
            else:
                raise StructureError('more than one atom satisfies the deletion condition')
        except StructureError:
            raise StructureError()
        
        
    def distance(self, atom1, atom2):
        """
        calculate the distance of two atom
        Note: don't consider the translational periodicity of crystal.
        
        Arguments:
            atom1/atom2: a array contained atomic information [element_symbol, x, y, z, coord_type (optional)]. i.e. ['Na', 0.5, 0.5, 0.5](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian'](Cartesian)
            
        Return:
            calculated distance between two atoms (float)
        """         
        # check the atom1 and atom2
        tmpa1, isExist1=self.isExist(atom1)
        tmpa2, isExist2=self.isExist(atom2)
        
        if not(isExist1) and not(isExist2):
            raise StructureError('non-exist atom1 and atom2!')
        elif not(isExist1):
            raise StructureError('non-exist atom1!')
        elif not(isExist2):
            raise StructureError('non-exist atom2!')
        
        # cartesian coordinate
        position1=direct2cartesian(self, extractingCoordinate(self, atom1))
        position2=direct2cartesian(self, extractingCoordinate(self, atom2))
        
        distance=np.linalg.norm(position2-position1)
        
        return distance
    
    def nearestNeighbors(self, atom, cutoff=3.0, **kwargs):
        """
        find out the nearest neighbors of a atom 

        Arguments:
            atom: a array contained atomic information [element_symbol, x, y, z, coord_type (optional)]. i.e. ['Na', 0.5, 0.5, 0.5](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian'](Cartesian)
            cutoff: the cutoff radius of searching neighbor
             
            kwargs:
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                tolerance (default=1e-3): distance tolerance between the position of the input and position of structure, unite is angstrom.

        Return:
            the nearest neighbors of a atom in [3,3,3] supercell by cartesian coordinate 
        """
        # step:
        # create supercell with the size [3,3,3] and make sure that the initial cell in the center of supercell cell.
        # according to the cutt-off radius, choose the neighbors of a atoms
        
        # check data
        # remove translation periodicity
        isNormalizingCoordinate=True
        if 'isNormalizingCoordinate' in kwargs:
            isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
            
        tmpa=self.exists(atom, isNormalizingCoordinate=isNormalizingCoordinate)

        size=[-1,2,-1,2,-1,2] # size of supercell
        a = size[1]-size[0]
        b = size[3]-size[2]
        c = size[5]-size[4]

        atom_set=self.atoms.all()
        satom_set=[] # atomic set of supercell
        for i in xrange(size[0], size[1]): # x
            for j in xrange(size[2], size[3]): # y
                for k in xrange(size[4], size[5]): # z
                    for atom in atom_set:
                        vector=i*self.lattice[0]+j*self.lattice[1]+k*self.lattice[2]
                        position = atom.direct2cartesian()
                        satom_set.append([atom.element.symbol, position[0]+vector[0], position[1]+vector[1], position[2]+vector[2]])

        #atom1_position=direct2cartesian(self,atom1[1:]) # the cartesian coordinate of the input atom1
        position0=tmpa.direct2cartesian() # position of given atom

        # the tolerance of coordination between the input atom and the atom in the structure, if it is smaller
        #                                            the value, two atoms are considered the same.
        tolerance = 1e-3 # default value
        if 'tolerance' in kwargs:
            tolerance=kwargs['tolerance']
            
        nn_set=[] # set of nearest neighbors
        for a in satom_set:
            distance=np.linalg.norm(a[1:]-position0)
            if (distance > tolerance) and (distance <= cutoff):
                position=cartesian2direct(self, a[1:])
                a=[a[0],position[0], position[1], position[2]]
                nn_set.append(a)
        
        return np.array(nn_set)
        
    def nearestNeighborsAll(self, cutoff=3.0, **kwargs):
        """
        list all nearest neighbors in a structure
        
        Arguments:
            cutoff: the cutoff radius of searching neighbor
             
            kwargs:
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                tolerance (default=1e-3): distance tolerance between the position of the input and position of structure, unite is angstrom.
        
        Return:
            
        """ 
        nna_set=[] # list of all nearest neighbors
        site_set=self.sites.all()
        for site in site_set:
            atom0=site.atoms.all()[0]
            nn_set=self.nearestNeighbors([atom0.element.symbol, atom0.x, atom0.y, atom0.z], cutoff).tolist()
            
            nna_set.append([[atom0.element.symbol, atom0.x, atom0.y, atom0.z],nn_set])
        
        return np.array(nna_set)
            
    def rdf(self, cutoff=10.0, delta=0.1):
        """
        radical distribution function
        
        
        """
        pass
    
    # add some method of pattern recognition?
    def islayered(self):
        """
        analysis if a structure is layered
    
        """
        pass
    