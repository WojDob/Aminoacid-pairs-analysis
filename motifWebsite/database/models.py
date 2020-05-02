from django.db import models

# Create your models here.

class Organism(models.Model):
    name = models.CharField(max_length=100)

class Motif(models.Model):
    motif = models.CharField(max_length=2)

class Protein(models.Model):
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE) # AZ
    proteinId = models.CharField(max_length=20)
    description = models.CharField(max_length=1000) # AZ (protein description from FASTA headers)
    length = models.PositiveIntegerField()
    sequence = models.TextField() # AZ

class Result(models.Model):
    motif = models.ForeignKey(Motif, on_delete=models.CASCADE)
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    _motif_positions = models.TextField() # Motif's positions in the protein sequence (e.g., "3,4,6")
    count = models.PositiveIntegerField()
    zscore = models.FloatField()
    ratio = models.FloatField()
    pvalue = models.FloatField()
    alignment = models.TextField() # AZ

    @property
    def motif_positions(self):
        # >>> results = Results.objects.get(id=1)
        # >>> print(results._motif_positions)
        # "3,4,6"
        # >>> print(results.motif_positions)
        # [3, 4, 6]
        return [int(pos) for pos in self._motif_positions.split(',')]
    