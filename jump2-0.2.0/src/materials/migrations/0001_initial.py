# -*- coding: utf-8 -*-
# Generated by Django 1.11.5 on 2018-01-09 20:16
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion
import utils.customField


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Atom',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ox', models.IntegerField(blank=True, null=True)),
                ('position', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('force', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('velocity', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('constraint', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('magmom', models.FloatField(blank=True, null=True)),
                ('charge', models.FloatField(blank=True, null=True)),
                ('volume', models.FloatField(blank=True, null=True)),
                ('occupancy', models.FloatField(default=1)),
            ],
            options={
                'db_table': 'atom',
                'default_related_name': 'atom_set',
            },
            bases=(models.Model, object),
        ),
        migrations.CreateModel(
            name='Case',
            fields=[
                ('name', models.CharField(max_length=255, primary_key=True, serialize=False)),
                ('calculated_parameters', utils.customField.DictField(blank=True, null=True)),
                ('energy', models.FloatField(blank=True, null=True)),
                ('energy_per_formula', models.FloatField(blank=True, null=True)),
                ('energy_per_atom', models.FloatField(blank=True, null=True)),
                ('pressure', models.FloatField(blank=True, null=True)),
                ('bandgap', models.FloatField(blank=True, null=True)),
                ('bandgap_img', models.ImageField(blank=True, null=True, upload_to=b'')),
                ('electron_mass', utils.customField.DictField(blank=True, null=True)),
                ('hole_mass', utils.customField.DictField(blank=True, null=True)),
            ],
            options={
                'db_table': 'case',
                'default_related_name': 'case_set',
            },
            bases=(models.Model, object),
        ),
        migrations.CreateModel(
            name='Composition',
            fields=[
                ('formula', models.CharField(max_length=255, primary_key=True, serialize=False)),
                ('generic', models.CharField(blank=True, max_length=255, null=True)),
                ('mass', models.FloatField(blank=True, null=True)),
            ],
            options={
                'db_table': 'composition',
                'default_related_name': 'composition_set',
            },
            bases=(models.Model, object),
        ),
        migrations.CreateModel(
            name='Element',
            fields=[
                ('symbol', models.CharField(max_length=4, primary_key=True, serialize=False)),
                ('z', models.IntegerField()),
                ('name', models.CharField(max_length=20)),
                ('group', models.IntegerField()),
                ('period', models.IntegerField()),
                ('mass', models.FloatField(null=True)),
                ('electronegativity', models.FloatField(null=True)),
            ],
            options={
                'ordering': ('z',),
                'db_table': 'element',
                'default_related_name': 'element_set',
            },
        ),
        migrations.CreateModel(
            name='MolAtom',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('position', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('magmom', models.FloatField(blank=True, null=True)),
                ('charge', models.FloatField(blank=True, null=True)),
                ('volume', models.FloatField(blank=True, null=True)),
                ('ox', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'molAtom',
                'default_related_name': 'atom_set',
            },
        ),
        migrations.CreateModel(
            name='MolComposition',
            fields=[
                ('formula', models.CharField(max_length=255, primary_key=True, serialize=False)),
                ('generic', models.CharField(blank=True, max_length=255, null=True)),
                ('mass', models.FloatField(null=True)),
            ],
            options={
                'db_table': 'molComposition',
                'default_related_name': 'composition_set',
            },
        ),
        migrations.CreateModel(
            name='MolElement',
            fields=[
                ('symbol', models.CharField(max_length=4, primary_key=True, serialize=False)),
                ('z', models.IntegerField()),
                ('name', models.CharField(max_length=20)),
                ('group', models.IntegerField()),
                ('period', models.IntegerField()),
                ('mass', models.FloatField(blank=True, null=True)),
                ('electronegativity', models.FloatField(blank=True, null=True)),
            ],
            options={
                'db_table': 'molElement',
                'default_related_name': 'element_set',
            },
        ),
        migrations.CreateModel(
            name='MolStructure',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('label', models.CharField(blank=True, max_length=80, null=True)),
                ('natoms', models.IntegerField(blank=True, null=True)),
                ('ntypes', models.IntegerField(blank=True, null=True)),
                ('volume', models.FloatField(blank=True, null=True)),
                ('volume_per_atom', models.FloatField(blank=True, null=True)),
                ('composition', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='structure_set', to='materials.MolComposition')),
                ('element_set', models.ManyToManyField(related_name='structure_set', to='materials.MolElement')),
            ],
            options={
                'db_table': 'molStructure',
                'default_related_name': 'structure_set',
            },
            bases=(models.Model, object),
        ),
        migrations.CreateModel(
            name='Operation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('operation', utils.customField.DictField(blank=True, null=True)),
            ],
            options={
                'db_table': 'operation',
                'default_related_name': 'operation_set',
            },
        ),
        migrations.CreateModel(
            name='Prototype',
            fields=[
                ('name', models.CharField(max_length=80, primary_key=True, serialize=False)),
                ('composition', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='materials.Composition')),
            ],
            options={
                'db_table': 'prototype',
            },
        ),
        migrations.CreateModel(
            name='Site',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('position', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('coordination_number', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'site',
                'default_related_name': 'site_set',
            },
            bases=(models.Model, object),
        ),
        migrations.CreateModel(
            name='Spacegroup',
            fields=[
                ('number', models.IntegerField(primary_key=True, serialize=False)),
                ('international', models.CharField(blank=True, max_length=20, null=True)),
                ('hm', models.IntegerField(blank=True, null=True)),
                ('hall', models.CharField(blank=True, max_length=20, null=True)),
                ('pearson', models.CharField(blank=True, max_length=20, null=True)),
                ('schoenflies', models.CharField(blank=True, max_length=20, null=True)),
                ('lattice_system', models.CharField(blank=True, max_length=20, null=True)),
                ('centerosymmetric', models.NullBooleanField()),
                ('operation_set', models.ManyToManyField(blank=True, related_name='spacegroup_set', to='materials.Operation')),
            ],
            options={
                'db_table': 'spacegroup',
                'default_related_name': 'spacegroup_set',
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('name', models.CharField(max_length=8, primary_key=True, serialize=False)),
                ('ox', models.IntegerField(blank=True, null=True)),
                ('element', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='species_set', to='materials.Element')),
            ],
            options={
                'db_table': 'species',
                'default_related_name': 'species_set',
            },
        ),
        migrations.CreateModel(
            name='Structure',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('label', models.CharField(blank=True, max_length=80, null=True)),
                ('comment', models.CharField(default='', max_length=80)),
                ('natoms', models.IntegerField(blank=True, null=True)),
                ('nsites', models.IntegerField(blank=True, null=True)),
                ('ntypes', models.IntegerField(blank=True, null=True)),
                ('multiple', models.IntegerField(blank=True, null=True)),
                ('lattice', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('volume', models.FloatField(blank=True, null=True)),
                ('volume_per_atom', models.FloatField(blank=True, null=True)),
                ('composition', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='structure_set', to='materials.Composition')),
                ('element_set', models.ManyToManyField(related_name='structure_set', to='materials.Element')),
                ('prototype', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='structure_set', to='materials.Prototype')),
                ('spacegroup', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='structure_set', to='materials.Spacegroup')),
                ('species_set', models.ManyToManyField(related_name='structure_set', to='materials.Species')),
            ],
            options={
                'db_table': 'structure',
                'default_related_name': 'structure_set',
            },
        ),
        migrations.CreateModel(
            name='WyckoffSite',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('symbol', models.CharField(max_length=1)),
                ('multiplicity', models.IntegerField(blank=True, null=True)),
                ('position', utils.customField.NumpyArrayField(blank=True, null=True)),
                ('spacegroup', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='wyckoffSite_set', to='materials.Spacegroup')),
            ],
            options={
                'db_table': 'wyckoffSite',
                'default_related_name': 'wyckoffSite_set',
            },
        ),
        migrations.AddField(
            model_name='site',
            name='structure',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='site_set', to='materials.Structure'),
        ),
        migrations.AddField(
            model_name='site',
            name='wyckoffSite',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='site_set', to='materials.WyckoffSite'),
        ),
        migrations.AddField(
            model_name='prototype',
            name='structure_of_prototype',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='+', to='materials.Structure'),
        ),
        migrations.AddField(
            model_name='molcomposition',
            name='element_set',
            field=models.ManyToManyField(related_name='composition_set', to='materials.MolElement'),
        ),
        migrations.AddField(
            model_name='molatom',
            name='element',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='atom_set', to='materials.MolElement'),
        ),
        migrations.AddField(
            model_name='molatom',
            name='structure',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='atom_set', to='materials.MolStructure'),
        ),
        migrations.AddField(
            model_name='composition',
            name='element_set',
            field=models.ManyToManyField(related_name='composition_set', to='materials.Element'),
        ),
        migrations.AddField(
            model_name='case',
            name='structure',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='case_set', to='materials.Structure'),
        ),
        migrations.AddField(
            model_name='atom',
            name='element',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='atom_set', to='materials.Element'),
        ),
        migrations.AddField(
            model_name='atom',
            name='site',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='atom_set', to='materials.Site'),
        ),
        migrations.AddField(
            model_name='atom',
            name='species',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='atom_set', to='materials.Species'),
        ),
        migrations.AddField(
            model_name='atom',
            name='structure',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='atom_set', to='materials.Structure'),
        ),
        migrations.AddField(
            model_name='atom',
            name='wyckoffSite',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='atom_set', to='materials.WyckoffSite'),
        ),
    ]