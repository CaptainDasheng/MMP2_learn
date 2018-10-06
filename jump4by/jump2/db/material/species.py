from __future__ import unicode_literals

from django.db import models

class Species(models.Model):
    # relationship
    element=models.ForeignKey('Element', blank=True, null=True)
    
    name=models.CharField(max_length=8, primary_key=True) # i.e. Na+
    ox=models.IntegerField(blank=True,null=True) # oxidation state
    
    class Meta:
        app_label='material'
        db_table='species'
        default_related_name='species'
    
    def __unicode__(self):
        return self.name
    