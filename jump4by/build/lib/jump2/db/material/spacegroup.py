from __future__ import unicode_literals

from django.db import models

# Create your models here.

class Translation(models.Model):
    # relationships
    #operations=models.ForeignKey('Operations', blank=True, null=True, related_name='translations')
    
    x=models.FloatField()
    y=models.FloatField()
    z=models.FloatField()
    
    class Meta:
        app_label='material'
        db_table='translation'
        default_related_name='translations'
        
    def __unicode__(self):
        return self.id
    
    
class Rotation(models.Model):
    # relationships
    #operations=models.ForeignKey('Operations', blank=True, null=True, related_name='rotations')
    
    a11=models.FloatField()
    a12=models.FloatField()
    a13=models.FloatField()
    a21=models.FloatField()
    a22=models.FloatField()
    a23=models.FloatField()
    a31=models.FloatField()
    a32=models.FloatField()
    a33=models.FloatField()
    
    class Meta:
        app_label='material'
        db_table='rotation'
        default_related_name='rotations'

    def __unicode__(self):
        return self.id
    
    
class Operation(models.Model):
    # relationships
    rotation=models.ForeignKey('Rotation', blank=True, null=True)
    translation=models.ForeignKey('Translation', blank=True, null=True)
    #spacegroups=models.ManyToManyField('Operation', null=True)
    
    class Meta:
        app_label='material'
        db_table='operation'
        default_related_name='operations'
        
    def __unicode__(self):
        return str(self.id)
    
    
class WyckoffSite(models.Model):
    # relationships
    spacegroup=models.ForeignKey('Spacegroup', blank=True, null=True)
    
    symbol=models.CharField(max_length=1)
    multiplicity=models.IntegerField(blank=True, null=True)
    x=models.FloatField(blank=True, null=True)
    y=models.FloatField(blank=True, null=True)
    z=models.FloatField(blank=True, null=True)
    
    class Meta:
        app_label='material'
        db_table='wyckoffSite'
        default_related_name='wyckoffSites'
        
    def __unicode__(self):
        return str(self.multiplicity)+self.symbol
            
            
class Spacegroup(models.Model):
    # relationships
    operations=models.ManyToManyField('Operation')
    
    number=models.IntegerField(primary_key=True)
    hm=models.CharField(max_length=20, blank=True, null=True)
    hall=models.CharField(max_length=20, blank=True, null=True)
    pearson=models.CharField(max_length=20, blank=True, null=True)
    schoenflies=models.CharField(max_length=20, blank=True, null=True)
    lattice_system=models.CharField(max_length=20, blank=True, null=True)
    centerosymmetric=models.BooleanField(default=False)
    
    class Meta:
        app_label='material'
        db_table='spacegroup'
        default_related_name='spacegroups'
    
    def __unicode__(self):
        return str(self.number)    