from __future__ import unicode_literals

import numpy as np
import spglib

from collections import OrderedDict

from django.db import models

from init_elements import InitializeElement

from prototype import Prototype

from structure import Structure
from composition import Composition
from element import Element
from species import Species
from atom import Atom, cartesian2direct, direct2cartesian

from spacegroup import Spacegroup

import sys, os
from __builtin__ import int
pro_dir=os.getcwd()
sys.path.append(pro_dir)
sys.path.append(pro_dir+'/../')
from utils.math import normalizingCoordinate, extractingCoordinate

class StructureFactoryError(Exception):
    pass


class StructureFactory(object):
    '''
    structure factory contain detecting self-attributes and fabricating new structures
    
    Arguments:
        structure: django's object of structure or structuring objects?
    
    '''
    def __init__(self, structure):
        self.__raw_structure=structure
        self.structure=structure.clone()
    
    def zoom(self, scale):
        """
        scale the lattice vector
        
        Arguments:
            scale: coefficient of zoom for lattice parameter
            
        Return:
            object of new structure
        """
        # note that: need to update the volume and volume_pa. check if other parameters also need to update.
        # method 1
        #self.structure.x1 *= scale
        #self.structure.x2 *= scale
        #self.structure.x3 *= scale
        #self.structure.y1 *= scale
        #self.structure.y2 *= scale
        #self.structure.y3 *= scale
        #self.structure.z1 *= scale
        #self.structure.z2 *= scale
        #self.structure.z3 *= scale
        # method 2
        lattice=self.structure.lattice
        self.structure.lattice=lattice*scale
        
        self.structure.volume=self.structure.calculatedVolume()
        
        return self.structure
        
    def addAtom(self, *atoms, **kwargs):
        """
        add some atoms to structure
        Note: considers translation periodicity of atomic coordinate.
        Arguments:
            atom: there are two format of the input for atom parameter:
                1. atom object
                2. a array contained atomic information [element_symbol, x, y, z, coord_type (optional)]. i.e. ['Na', 0.5, 0.5, 0.5](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian'](Cartesian)
         
        kwargs:
                isConstraint:True/False (default=False).
                    
        Return:
            object of new structure
        """
        # check length of atoms
        if len(atoms) == 0:
            raise StructureFactoryError("atoms can't be []")
        
        for atom in atoms:
            
            # check type of atom
            if isinstance(atom, Atom):
                if atom.structure == self.__raw_structure:
                    #raise StructureFactoryError("atom already exists in the structure")
                    self.structure.pop()
                    return self.__raw_structure
                    break
                atom=[atom.element.symbol, atom.x, atom.y, atom.z]
        
            # check whether existing this atom on this position
            if atom[1:] in self.structure: # exists this site
                self.structure.pop()
                return self.__raw_structure
                break
            else:
                # check element
                element_symbol=atom[0] # need to check valid
                if Element.objects.filter(symbol=element_symbol).count() == 0: # nonexistent element in element table
                    InitializeElement(element_symbol)
                if not(element_symbol in self.structure): # nonexistent element in the structure
                    self.structure.elements.add(Element.objects.get(symbol=element_symbol))
                
                # position
                position=normalizingCoordinate(extractingCoordinate(self.structure, atom))
        
                # constraint
                constrained=False # whether constrain atom
                if 'isConstraint' in kwargs:
                    constrained=kwargs['isConstraint']
            
                atom=Atom.objects.create(structure=self.structure,
                                         element=Element.objects.get(symbol=element_symbol),
                                         x=position[0],
                                         y=position[1],
                                         z=position[2])

        poscar=self.structure.formatting(isConstraint=constrained)
        self.structure.pop() # delete old object
        self.structure=Structure().append(poscar)
        
        return self.structure

    def deleteAtom(self, *atoms, **kwargs):
        """
        delete some atoms
        
        Arguments:
            atom: there are two format of the input for atom parameter:
                1. atom object
                2. a array contained atomic information [element_symbol, x, y, z, coord_type (optional)]. i.e. ['Na', 0.5, 0.5, 0.5](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian'](Cartesian)
              
        
            kwargs:
                tolerance (default=1e-3): distance tolerance between the position of the input and position of structure, unite is angstrom.
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                isConstraint:True/False (default=False).
        
        Return:
            object of new structure and status of operating (True/False).
        """
        # check data
        # check length of atoms
        if len(atoms) == 0:
            raise StructureFactoryError("atoms can't be []")
        
        for atom in atoms:        

            # check type of atom
            if isinstance(atom, Atom):
                if atom.structure != self.__raw_structure:
                    raise StructureFactoryError("atom doesn't belong the structure")
                atom=[atom.element.symbol, atom.x, atom.y, atom.z]
            
            # tolerance
            tolerance=1e-3
            if 'tolerance' in kwargs:
                tolerance=kwargs['tolerance']
            
            # remove translation periodicity
            isNormalizingCoordinate=True
            if 'isNormalizingCoordinate' in kwargs:
                isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
            # constraint
            constrained=False # whether constrain atom
            if 'isConstraint' in kwargs:
                constrained=kwargs['isConstraint']
            
            tmpa=self.structure.exists(atom, isNormalizingCoordinate=isNormalizingCoordinate, isConstraint=constrained)    
            
            if tmpa == None:
                self.structure.pop()
                return self.__raw_structure
                break
            else:
                tmpa.delete()
                
        poscar=self.structure.formatting(isConstraint=constrained)
        self.structure.pop() # delete old object
        self.structure=Structure().append(poscar)
            
        return self.structure    
    
    def substituteAtom(self, *atoms, **kwargs):
        """
        change type (element or species?) of a atom
        
        Arguments:
            atoms: there are two format of the input for atom parameter:
                1. [atom_object, new_element_symbol]
                2. [element_symbol, x, y, z, coord_type (optional), new_element_symbol]. i.e. ['Na', 0.5, 0.5, 0.5, 'K'](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian', 'K'](Cartesian)
            Note: the format of atoms is different with that of atoms in addAtoms and deleteAtoms methods.
            
            #atom: there are two format of the input for atom parameter:
            #    1. atom object
            #    2. a array contained atomic information [element_symbol, x, y, z, coord_type (optional)]. i.e. ['Na', 0.5, 0.5, 0.5](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian'](Cartesian)
            
            kwargs:
                tolerance (default=1e-3): distance tolerance between the position of the input and position of structure, unite is angstrom.  
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                isConstraint:True/False (default=False).
                
        Return:
            object of new structure                
        """
        # check data
        # check length of atoms
        if len(atoms) == 0:
            raise StructureFactoryError("atoms can't be []")
        
        for atom in atoms:        

            if len(atom) ==2:
                element_symbol=atom[1]
                atom=atom[0]
            else:
                element_symbol=atom[-1]
                atom=atom[:-1]
                
            # check type of atom
            if isinstance(atom, Atom):
                if atom.structure != self.__raw_structure:
                    raise StructureFactoryError("atom doesn't belong the structure")
                atom=[atom.element.symbol, atom.x, atom.y, atom.z]
        
            # tolerance
            tolerance=1e-3
            if 'tolerance' in kwargs:
                tolerance=kwargs['tolerance']
        
            # remove translation periodicity
            isNormalizingCoordinate=True
            if 'isNormalizingCoordinate' in kwargs:
                isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
            # constraint
            constrained=False # whether constrain atom
            if 'isConstraint' in kwargs:
                constrained=kwargs['isConstraint']
            
            tmpa=self.structure.exists(atom, isNormalizingCoordinate=isNormalizingCoordinate, isConstraint=constrained)    
            
            if tmpa == None:
                self.structure.pop()
                return self.__raw_structure
                break
            else:
                # check element_symbol
                if Element.objects.filter(symbol=element_symbol).count() == 0: # nonexistent element in element table
                    InitializeElement(element_symbol)
                    
                # add element to structure
                self.structure.elements.add(Element.objects.get(symbol=element_symbol))
                self.structure.save()                    
                    
                # substitute
                # method 1 (1.create new atom, 2. delete old atom)
                atom=Atom.objects.create(structure=self.structure,
                                         element=Element.objects.get(symbol=element_symbol),
                                         x=tmpa.position[0],
                                         y=tmpa.position[1],
                                         z=tmpa.position[2])
                tmpa.delete()
                atom.save()
                
                # method 2 (modify atom.element)
                # atom_set[index].element=Element.objects.get(symbol=element_symbol)
                # atom_set[index].save()
                
        poscar=self.structure.formatting(isConstraint=constrained)
        self.structure.pop() # delete old object
        self.structure=Structure().append(poscar)
            
        return self.structure        
        
    def vacuum(self, direction, **kwargs):
        """
        add vacuum along a direction
        
        arguments:
            direction: direction vector to add the vacuum along lattice vector(a/b/c) (unit: Angstroms). i.e. [2, 0, 0]
            
            kwargs:
                isCenter: centroid of all atoms in structure (True/False). The default value is True.
                distance: moving distance for all atoms. i.e. 2. Note that distance <= modulus of direction vector
                isConstraint:True/False (default=False).
        Return:
            object of new structure
        """
        # check direction, ensure to exist two zero and one non-zero.
        if len(direction) != 3:
            raise StructureFactoryError('invalid dimension of direction!')
        
        # method 1
        #dcounter={x:direction.count(x) for x in direction}
        #if dcounter[0] > 2:
        # method 2
        dsum=np.sum(direction) # sum of direction
        dabsum=np.sum(np.absolute(direction)) # sum of absolute value of direction
        if (dsum != dabsum) or dabsum == 0:
            raise StructureFactoryError('invalid direction!')
        
        isCenter=True # default
        if 'isCenter' in kwargs:
            isCenter=kwargs['isCenter']
            if 'distance' in kwargs:
                raise StructureFactoryError('conflicts between isCenter and distance parameters! only exits one.')
        distance=0
        if 'distance' in kwargs:
            distance=kwargs['distance']
        
        # constraint
        constrained=False # whether constrain atom
        if 'isConstraint' in kwargs:
            constrained=kwargs['isConstraint']
        
        # transform coordinate of direction to cartesian.
        lattice=self.structure.lattice
        lattice_parameters=self.structure.lattice_parameters
        dtmp=[]
        for i in xrange(0, len(direction)):
            if direction[i] != 0: # vacuum direction
                dtmp=np.array(lattice[i]*direction[i]/lattice_parameters[i]) # new translation direction

                # move atoms
                atom_set=self.structure.atoms.all()
                atom_position_set={} # temporary save the new position of atom. {atom:new_position(cartesian)}
                atom_position_set2={} # temporary save the new position of atom which its component equal 0 along vacuum direction.
                for j in xrange(0, atom_set.count()):
                    atom=atom_set[j]
                    atom_position_set[atom]=atom.direct2cartesian()
                    
                    if atom.position[i] == 0:
                        
                        position=atom.position
                        position[i]=1
                        atmp=Atom.objects.create(structure=self.structure,
                                                 element=atom.element,
                                                 x=position[0],
                                                 y=position[1],
                                                 z=position[2])
                        atom_position_set2[atmp]=atom.direct2cartesian()+lattice[i]
   
                # change lattice
                scale=(self.structure.lattice_parameters[i]+direction[i])/(self.structure.lattice_parameters[i])
                if i == 0: # a
                    self.structure.x1=(self.structure.x1)*scale
                    self.structure.y1=(self.structure.y1)*scale
                    self.structure.z1=(self.structure.z1)*scale
                elif i == 1: # b
                    self.structure.x2=(self.structure.x2)*scale
                    self.structure.y2=(self.structure.y2)*scale
                    self.structure.z2=(self.structure.z2)*scale
                elif i == 2: # c
                    self.structure.x3=(self.structure.x3)*scale
                    self.structure.y3=(self.structure.y3)*scale
                    self.structure.z3=(self.structure.z3)*scale
                self.structure.save()

                # change atomic position
                # distance
                if 'distance' in kwargs:
                    scale=distance/(self.structure.lattice_parameters[i]+direction[i])
                    for atom, position in atom_position_set.iteritems():
                        atom.position=np.array(cartesian2direct(self.structure, position+dtmp*scale))
                    for atom, position in atom_position_set2.iteritems():
                        atom.position=np.array(cartesian2direct(self.structure, position+dtmp*scale))
                # isCenter                
                elif isCenter:
                    for atom, position in atom_position_set.iteritems():
                        atom.position=np.array(cartesian2direct(self.structure, position+dtmp/2))
                    for atom, position in atom_position_set2.iteritems():
                        atom.position=np.array(cartesian2direct(self.structure, position+dtmp/2))
                else: # isCenter=False, atom don't move
                    for atom, position in atom_position_set.iteritems():
                        atom.position=np.array(cartesian2direct(self.structure, position))
                    for atom, position in atom_position_set2.iteritems():
                        atom.position=np.array(cartesian2direct(self.structure, position))
        
        poscar=self.structure.formatting(isConstraint=constrained)
        self.structure.pop() # delete old object
        self.structure=Structure().append(poscar)
                
        return self.structure
    
    def magnetismOrder(self, elements):
        """
        At present, only consider FM configuration. Other magnetic configuration need to set the atomic magnetism one by one.
        
        Arguments
            elements: dict of element's symbol and its magnetic moment.
                {'Fe': 5,
                 'Cr': 3,
                     ...}
        """
        for k in elements:
            if self.structure.exists(k) == None:
                raise StructureFactoryError('nonexistent element')
        
        for atom in self.structure.atoms.all():
            if atom.element.symbol in elements:
                atom.magmom=elements[atom.element.symbol]
            else:
                atom.magmom=0
            atom.save()
                
        return self.structure
    
    def constraint(self, *atoms, **kwargs):
        """
        selected dynamics
        
        Arguments:
            atom: there are two format of the input for atom parameter:
                1. [atom_object, constraint]
                2. a array contained atomic information [element_symbol, x, y, z, coord_type (optional), constraint]. i.e. ['Na', 0.5, 0.5, 0.5, [True, True, False]](Direc) or ['Na', 1.3, 3.5, 2.0, 'Cartesian', [True, True, False]](Cartesian)
                
                constraint: constraint of atomic position (boolean value). i.e. [True, True, False]
            
            kwargs:
                tolerance: distance tolerance between the position of the input and position of structure, unite is angstrom.  
                isNormalizingCoordinate (default=True): consider the translation periodicity for input of atomic coordinate when calculating the distance, ensure the value is between 0 and 1. i.e. 2.3 -> 0.3.
                isConstraint:True/False (default=False).
        
        Return:
            object of new structure and status of operating (True/False).       
        """
        # check data
        # check length of atoms
        if len(atoms) == 0:
            raise StructureFactoryError("atoms can't be []")
        
        for atom in atoms:        

            if len(atom) ==2:
                constraint=atom[1]
                atom=atom[0]
            else:
                constraint=atom[-1]
                atom=atom[:-1]
                
            # check type of atom
            if isinstance(atom, Atom):
                if atom.structure != self.__raw_structure:
                    raise StructureFactoryError("atom doesn't belong the structure")
                tom=[atom.element.symbol, atom.x, atom.y, atom.z]
        
            # tolerance
            tolerance=1e-3
            if 'tolerance' in kwargs:
                tolerance=kwargs['tolerance']
            
            # remove translation periodicity
            isNormalizingCoordinate=True
            if 'isNormalizingCoordinate' in kwargs:
                isNormalizingCoordinate=kwargs['isNormalizingCoordinate']
            # constraint
            constrained=False # whether constrain atom
            if 'isConstraint' in kwargs:
                constrained=kwargs['isConstraint']
            
            tmpa=self.structure.exists(atom, isNormalizingCoordinate=isNormalizingCoordinate, isConstraint=constrained)
        
            if tmpa == None:
                self.structure.pop()
                return self.__raw_structure
                break
            else:
                tmpa.cx = constraint[0]
                tmpa.cy = constraint[1]
                tmpa.cz = constraint[2]
                tmpa.save()

        poscar = self.structure.formatting(isConstraint=True)
        self.structure.pop()  # delete old object
        self.structure = Structure().append(poscar)

        return self.structure

    def redefine(self, operatorMatrix):
        """
        redefine lattice cell. C'=C x M 
        
        Arguments:
            operatorMatrix: operator matrix (i.e. M).
                M=[[0, 1, 1],
                   [1, 0, 1],
                   [1, 1, 0]]
            Note: the component of M should be integer. And the volume of M is an integer greater than 0.
        
        Return:
            object of structure and status of operating
        """
        # check operatorMatrix
        tmpm=np.array(operatorMatrix)
        if not((len(tmpm.shape) == 2) and ((tmpm[0] == 3) and (tmpm[1] == 3))):
            raise StructureFactoryError('invalid dimension of operatorMatrix (3x3)')
        
        for component in np.nditer(tmpm):
            if not(isinstance(component, int)):
                raise StructureFactoryError("component of operatorMatrix isn't an integer")
        volume=np.linalg.det(tmpm)
        if not(isinstance(volume, int)):
            raise StructureFactoryError("volume of operatorMatrix isn't an integer")
        
        
        
    def __toPOSCAR(self, cell):
        """
        convert cell to poscar   
        """
        lattice=cell[0]

        # elements, numbers and position
        atom_set=OrderedDict() # elements and its positions
        
        for i in xrange(0, len(cell[2])):
            e=Element.objects.filter(z=cell[2][i])[0].symbol
            if e in atom_set:
                tmp=atom_set[e]
                tmp.append(cell[1][i].tolist())
                atom_set[e]=tmp
            else:
                atom_set[e]=[cell[1][i].tolist()]
        
        elements=[]
        numbers=[]
        positions=[]
        for k in atom_set.keys():
            elements.append(k)
            numbers.append(len(atom_set[k]))
            if positions == []:
                positions=atom_set[k]
            else:
                positions=np.vstack((positions, atom_set[k]))

        elements = np.array(elements)
        numbers = np.array(numbers)

        # poscar=(comment,lattice,elements,numbers,type,positions,constraints)
        poscar = {'lattice': lattice,
                  'elements': elements,
                  'numbers': numbers,
                  'type': 'Direct',
                  'positions': positions}
        
        return poscar
    
    def refine(self, symprec=1e-5, angle_tolerance=-1.0):
        """
        refine structure which can change the cell's shape
        
        Arguments:
            symprec (symmetry tolerance): distance tolerance in Cartesian coordinates to find crystal symmetry
            angle_tolerance: An experimental argument that controls angle tolerance between basis vectors.
                Normally it is not recommended to use this argument.
    
        """
        cell=self.structure.formatting('cell')
        cell = (cell['lattice'], cell['positions'], cell['numbers'])
        newcell=spglib.refine_cell(cell, symprec=symprec, angle_tolerance=angle_tolerance)
        
        if newcell == None:
            raise StructureFactoryError('the search is filed')
        else:
            poscar=self.__toPOSCAR(newcell)
            
            self.structure.pop()  # delete old object
            self.structure = Structure().append(poscar)
            return self.structure
    
    def primitive(self, symprec=1e-5):
        """
        primitive structure
        
        Arugments:
            symprec (symmetry tolerance): distance tolerance in Cartesian coordinates to find crystal symmetry 
        """
        cell=self.structure.formatting('cell')
        cell = (cell['lattice'], cell['positions'], cell['numbers'])
        newcell=spglib.find_primitive(cell, symprec=symprec)
        
        if newcell == None:
            raise StructureFactoryError('the search is filed')
        else:
            poscar=self.__toPOSCAR(newcell)
            
            self.structure.pop()  # delete old object
            self.structure = Structure().append(poscar)
            return self.structure
    
    def conventional(self, symprec=1e-5):
        """
        conventional structure
        
        Arugments:
            symprec (symmetry tolerance): distance tolerance in Cartesian coordinates to find crystal symmetry 
        """
        cell=self.structure.formatting('cell')
        cell = (cell['lattice'], cell['positions'], cell['numbers'])
        
        # method 1
        lattice, scaled_positions, numbers=spglib.standardize_cell(cell, symprec=symprec)
        newcell=(lattice, scaled_positions, numbers)
        # method 2
        #newcell=spglib.standardize_cell(cell, symprec=symprec)
        
        if newcell == None:
            raise StructureFactoryError('The search is failed')
        else:
            poscar=self.__toPOSCAR(newcell)
            
            self.structure.pop()  # delete old object
            self.structure = Structure().append(poscar)
            return self.structure
    
    def niggliReduce(self, eps=1e-5):
        """
        Niggli reduction
        
        Arguments:
            eps: tolerance parameter, but unlike symprec the unit is not a length.
                This is used to check if difference of norms of two basis vectors 
                is close to zero or not and if two basis vectors are orthogonal
                by the value of dot product being close to zero or not. The detail
                is shown at https://atztogo.github.io/niggli/.
        """
        cell=self.structure.formatting('cell')
        lattice=cell['lattice']
        niggli_lattice=spglib.niggli_reduce(lattice, eps=eps)
        
        return niggli_lattice
    
    def delaunayReduce(self, eps=1e-5):
        """
        Delaunay reduction
        
        Arguments:
            eps: tholerance parameter, see niggliReduce.
        """
        cell=self.structure.formatting('cell')
        lattice=cell['lattice']
        delaunay_lattice=spglib.delaunay_reduce(lattice, eps=eps)
    
        return delaunay_lattice
    
    def supercell(self, dim):
        """
        create supercell
        
        Arguments:
            dim: size of supercell. i.e. [2,2,2]
            
        Return:
            object of supercell structure.
        """
        # step:
        # 1. move atom
        # 2. change the lattice vector
        # 3. call formatting
        # 4. pop old
        # 5. append new
        
        # check dim
        if (sum(1 for x in dim if x <= 0) >= 1) or (sum(1 for x in dim if not isinstance(x, int)) >= 1):
            raise StructureFactoryError('invalid dim')
        
        atom_set=self.structure.atoms.all()
        for i in xrange(0, dim[0]): # x
            for j in xrange(0, dim[1]): # y
                for k in xrange(0, dim[2]): # z
                    for atom in atom_set:
                        Atom.objects.create(structure=self.structure,
                                            element=atom.element,
                                            x=(atom.x+i)/dim[0],
                                            y=(atom.y+j)/dim[1],
                                            z=(atom.z+k)/dim[2])
        # delete old atoms
        for atom in atom_set:
            atom.delete()

        lattice=self.structure.lattice
        for i in xrange(0, lattice.shape[0]):
            lattice[i]=lattice[i]*dim[i]
                   
        self.structure.lattice=lattice
        self.structure.save()
        
        #poscar=self.structure.formatting(isConstraint=constrained)
        poscar=self.structure.formatting(isConstraint=False)
        self.structure.pop() # delete old object
        self.structure=Structure().append(poscar)

        return self.structure
    
    def alloy(self):
        """
        create alloy
        
        
        """
        pass
    
    def surface(self):
        """
        
        """
        
    def adsorption(self):
        """
        adsorption on a surface
        
        
        """
        pass
    
    def rotation(self, axis, angle, angle_type="Radian", **kwargs):
        """
        rotate selected atoms set
        
        Arguments:
            axis: rotating axis [x, y, z].
            angle: rotating angle.
            angle_type: type of angle (Radian(default)/Degree).

            kwargs:
                isConstraint:True/False (default=False).

        Return:
            object of new structure
        """
        axis=axis/np.linalg.norm(axis)
        rx=axis[0]
        ry=axis[1]
        rz=axis[2]

        type='Radian'
        if angle_type.lower().startswith('d'):
            type='Degree'
            angle=np.deg2rad(angle)
        elif angle_type.lower().startswith('r'):
            type='Radian'
        else:
            raise StructureFactoryError('unknown value in angle_type(Radian/Degree)!')

        # constraint
        constrained=False  # whether constrain atom
        if 'isConstraint' in kwargs:
            constrained=kwargs['isConstraint']

        sina=np.sin(angle)
        cosa=np.cos(angle)

        #set the transition matrix
        r_matrix=np.matrix([[cosa+(1-cosa)*np.square(rx), (1-cosa)*rx*ry-rz*sina, (1-cosa)*rx*rz+ry*sina],
                            [(1-cosa)*rx*ry+rz*sina, cosa+(1-cosa)*np.square(ry), (1-cosa)*ry*rz-rx*sina],
                            [(1-cosa)*rx*rz-ry*sina, (1-cosa)*ry*rz+rx*sina, cosa+(1-cosa)*np.square(rz)]])

        #rotate the atoms
        for atom in self.structure.atoms.all():
            position=atom.position
            position=(position*r_matrix).tolist()[0]
            atom.x=position[0]
            atom.y=position[1]
            atom.z=position[2]
            atom.save()

        poscar=self.structure.formatting(isConstraint=constrained)
        self.structure.pop()  # delete old object
        self.structure=Structure().append(poscar)

        return self.structure
    
    def translation(self, atom_set, destination):
        """
        translation selected atoms set
        
        Arguments:
            atom_set: collection of operated atoms
            
                [[element_symbol1, x1, y1, z1, coord_type1 (optional; default=Direct)],
                 [element_symbol2, x2, y2, z2, coord_type2 (optional; default=Direct)],
                                                             ...]]
                
            ##vector: distance of translation(unit: Angstroms). i.e. [5.0, 0.0, 0.0]
            destination: position of moving destination.
                for crystal:
                    [x, y, z, coord_type (optional; default=Direct)]
                for molecule:
                    [x, y, z] # only Cartesian
        """
        if isinstance(self.structure, Structure): # crystal
            if (len(destination) == 3) or ((len(destination) == 4) and destination.lower().startswith('d')): # direct
                dposition=direct2cartesian(self.structure, destination) # position of destination
                
                
        elif isinstance(self.structure, MolStructure): # molecule
            pass

    def images4NEB(self, struc2, nstruc, **kwargs):
        """
        image structures for neb calculation
        
        Note that: beware whether the atomic orders between two structures are consistent.
        
        Arguments:
            struc2: the last structure
            nstruc: number of images(including the two endpoint structures)
            
            kwargs:
                isConstraint (default=True):True/False
                isTranslation (default=True): consider the translation periodicity when calculating the distance
        Return:
            a list of all images (structures).
        """
        if not ((self.structure.lattice.tolist() == struc2.lattice.tolist()) or (self.structure.natoms == struc2.natoms)):
            raise StructureFactoryError('the two structures are mismatched, check lattice vectors or number of atoms')
        
        isConstraint=True # default
        if 'isConstraint' in kwargs:
            isConstraint=kwargs['isConstraint']
            
        isTranslation=True
        if 'isTranslation' in kwargs:
            isTranslation=kwargs['isTranslation']
        
        atom_set1=self.structure.atoms.all() # atoms of the start structure
        atom_set2=struc2.atoms.all() # atoms of the end structure
        struc_set=[] # output structure list
        dvector_set=[] # atomic displacement (delta vector) for every images from the start structure to the end structure
        
        # vectors per step
        for i in xrange(0, len(atom_set1)):
            vector=[(atom_set2[i].x-atom_set1[i].x),
                     (atom_set2[i].y-atom_set1[i].y),
                     (atom_set2[i].z-atom_set1[i].z)]
            # translation periodicity
            if isTranslation:
                for j in xrange(0, len(vector)):
                    if vector[j] > 0.5:
                        vector[j] -= 1
                    elif vector[j] < -0.5:
                        vector[j] += 1
            
            dvector=np.array(vector)/(nstruc-1)
            dvector_set.append(dvector)

        # append the start structure
        struc_set.append(self.structure)

        # append the following structures
        tmps=self.structure.clone() # temporary structure
        for i in xrange(1, nstruc-1):
            atom_index=0
            for atom in tmps.atoms.all():
                atom.x += dvector_set[atom_index][0]
                atom.y += dvector_set[atom_index][1]
                atom.z += dvector_set[atom_index][2]
                atom.save()
                atom_index += 1
            poscar=tmps.formatting(isConstraint=isConstraint)
            struc_set.append(Structure().append(poscar))
        tmps.pop() # delete temporary structure
        
        # append the end structure
        struc_set.append(struc2)

        return struc_set

    
    def defect(self):
        """
        """
        pass
    
    def interface(self):
        """
        """
        pass
    
