from __future__ import unicode_literals
import numpy as np
from django.db import models

    
def cartesian2direct(structure, position):
    """
    cartesian to direct
    
    arguments:
        structure: object of structure
        position: a arrary of this atomic position (cartesian) (i.e. [x, y, z]).
        
    return:
        a arrary of this atomic position (direct) (i.e. [x, y, z]).
    """
    lattice=structure.lattice

    inv=np.linalg.inv(lattice).T
    direct=np.dot(inv,position)

    return direct


def direct2cartesian(structure, position):
    """
    direct to cartesian
    
    arguments:
        structure: object of structure
        position: a arrary of this atomic position (direct) (i.e. [x, y, z]).
        
    return:
        a arrary of this atomic position (cartesian) (i.e. [x, y, z]).
    """
    lattice=structure.lattice
    # method 1
    #cartesian=0
    #for i in xrange(0, 3): # dirction
    #    cartesian += position[i]*lattice[i]
      
    # method 2
    cartesian=np.dot(position, lattice)
    
    return np.array(cartesian)

class AtomError(Exception):
    pass


class Atom(models.Model):
    
    class Meta:
        unique_together = (('element', 'x', 'y', 'z'),)
    
    # relationships
    structure=models.ForeignKey('Structure', null=True)
    site=models.ForeignKey('Site', null=True)
    
    # species
    element=models.ForeignKey('Element', null=True)
    species=models.ForeignKey('Species', null=True)
    ox=models.IntegerField(default=None, blank=True, null=True) # oxidation state
    
    # position
    x=models.FloatField()
    y=models.FloatField()
    z=models.FloatField()
    
    # forces
    fx=models.FloatField(blank=True, null=True)
    fy=models.FloatField(blank=True, null=True)
    fz=models.FloatField(blank=True, null=True)
    
    # constraint
    cx=models.BooleanField(default=True)
    cy=models.BooleanField(default=True)
    cz=models.BooleanField(default=True)
    
    # properties
    magmom=models.FloatField(blank=True, null=True)
    charge=models.FloatField(blank=True, null=True) # effective valence
    volume=models.FloatField(blank=True, null=True)
    
    # symmetry
    occupancy=models.FloatField(default=1)
    wyckoff=models.ForeignKey('WyckoffSite', blank=True, null=True)
    
    class Meta:
        app_label='material'
        db_table='atom'
        default_related_name='atoms'
    
    def __unicode__(self):
        return '%s %f %f %f' %(self.element, self.x, self.y, self.z)
    
    # override django's method
    def delete(self, *args, **kwargs):
        # delete atom
        super(Atom, self).delete(*args, **kwargs) # Call the 'real' delete()
        
        # check number of atoms corresponding to its element. Need to delete the element when doesn't exist atom. 
        if self.element.atoms.filter(structure=self.structure).count() == 0:
            self.element.delete()
    
    @property
    def position(self):
        self._position=np.array([self.x, self.y, self.z])
        return self._position
    
    @position.setter
    def position(self, position):
        if len(position) != 3:
            raise AtomError('invalid position!')
        self.x=position[0]
        self.y=position[1]
        self.z=position[2]
        self.save()
        
    @property    
    def force(self):
        self._force=np.array([self.fx, self.fy, self.fz])
        return self._force
    
    @property
    def constraint(self):
        self._constraint=np.array([self.cx, self.cy, self.cz])
        return self._constraint
    
    @constraint.setter
    def setConstrain(self, constraint):
        self.cx=constraint[0]
        self.cy=constraint[1]
        self.cz=constraint[2]
        
    @property
    def coordiantion_number(self):
        return self.site.coordiantion_number

    # Note: in whole jump2, coordinate type of atom is direct.
    #def cartesian2direct(self):
    #    """
    #    cartesian to direct
    #    """
    #    return cartesian2direct(self.structure, self.position)
        
    def direct2cartesian(self):
        """
        direct to cartesian
        
        return:
            a arrary of this atomic position (cartesian) (i.e. [x, y, z]).
        """
        return direct2cartesian(self.structure, self.position)
    
    
    
    