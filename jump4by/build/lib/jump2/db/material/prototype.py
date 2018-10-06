from __future__ import unicode_literals

from django.db import models

class Prototype(models.Model):
    # relationships
    composition=models.ForeignKey('Composition', blank=True, null=True)
    structures=models.ForeignKey('Structure', blank=True, null=True, related_name='+')
    
    name=models.CharField(max_length=80, primary_key=True)
    
    class Meta:
        app_label='material'
        db_table='prototype'
        default_related_name='prototypes'