from __future__ import unicode_literals

from django.db import models

class MolComposition(models.Model):
    """
    Composition
    
    Relationships:
        Element
        Structure
        Prototype
    
    Attrubutes:
        formula: normalized composition. i.e. PbTiO3
        generic: genericized composition. i.e. ABO3
        mass : mass per formula
    """
    
    # relationship
    elements=models.ManyToManyField('MolElement')
    
    formula=models.CharField(primary_key=True, max_length=255)
    generic=models.CharField(max_length=255, blank=True, null=True) # i.e. ABO3
    mass=models.FloatField(null=True)
    
    class Meta:
        app_label='material'
        db_table='molComposition'
        default_related_name='compositions'
        
    def __unicode__(self):
        return self.formula
    
    
    