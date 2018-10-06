from __future__ import unicode_literals

from django.db import models

class Site(models.Model):
    # relationships
    structure=models.ForeignKey('Structure', blank=True, null=True)
    wyckoff=models.ForeignKey('WyckoffSite', blank=True, null=True)
    
    # positions
    x=models.FloatField()
    y=models.FloatField()
    z=models.FloatField()
    
    # coordination number
    coordiantion_number=models.IntegerField(blank=True, null=True)
    
    class Meta:
        app_label='material'
        db_table='site'
        default_related_name='sites'
    
    def __unicode__(self):
        #return '%f %f %f' %(self.x, self.y, self.z) #, self.wyckoff.)
        return '%s (%s) %f %f %f' %(self.atoms.filter(structure=self.structure)[0].element.symbol,
                                    str(self.wyckoff.multiplicity)+self.wyckoff.symbol,
                                    self.x, self.y, self.z)