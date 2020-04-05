from django.db import models

# Create your models here.

class Organism(models.Model):
	name = models.CharField(max_length=100)

class Motif(models.Model):
	motif = models.CharField(max_length=2)

class Protein(models.Model):
	origin = models.ForeignKey(Organism, on_delete=models.CASCADE)
	proteinId = models.CharField(max_length=20)
	length = models.PositiveIntegerField()

class Result(models.Model):
	motif = models.ForeignKey(Motif, on_delete=models.CASCADE)
	protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
	count = models.PositiveIntegerField()
	zscore = models.FloatField()
	ratio = models.FloatField()
	pvalue = models.FloatField()
	alignment = models.CharField(max_length=1000000)

