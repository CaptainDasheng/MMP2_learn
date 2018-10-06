from __future__ import unicode_literals

from django.db import models

class Element(models.Model):
    symbol=models.CharField(primary_key=True, max_length=4)
    z=models.IntegerField()
    name=models.CharField(max_length=20)
    
    # Periodic table
    group=models.IntegerField()
    period=models.IntegerField()
    
    # Physical characteristics
    mass=models.FloatField(null=True)
    
    # chemical characteristics
    electronegativity=models.FloatField(null=True)
    
    class Meta:
        app_label='material'
        db_table='element'
        default_related_name='elements'
            
    def __unicode__(self):
        return self.symbol
    