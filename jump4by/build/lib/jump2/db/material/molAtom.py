from __future__ import unicode_literals
import numpy as np
from django.db import models


class MolAtomError(Exception):
    pass


class MolAtom(models.Model):
    
    class Meta:
        unique_together = (('element', 'x', 'y', 'z'),)
    
    # relationships
    structure=models.ForeignKey('MolStructure', null=True)
    
    # species
    element=models.ForeignKey('MolElement', null=True)
    ox=models.IntegerField(default=None, blank=True, null=True) # oxidation state
    
    # position
    x=models.FloatField()
    y=models.FloatField()
    z=models.FloatField()
    
    # properties
    magmom=models.FloatField(blank=True, null=True)
    charge=models.FloatField(blank=True, null=True) # effective valence
    volume=models.FloatField(blank=True, null=True)
    
    
    class Meta:
        app_label='material'
        db_table='molAtom'
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
            raise MolAtomError('invalid position!')
        self.x=position[0]
        self.y=position[1]
        self.z=position[2]
        self.save()

        