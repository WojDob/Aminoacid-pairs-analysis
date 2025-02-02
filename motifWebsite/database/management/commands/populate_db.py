from django.core.management.base import BaseCommand
from database.models import *
from Bio import SeqIO
import json


class Command(BaseCommand):
	args = 'idontknow'
	help = 'yes please'

	def _create_organisms(self):
		athaliana = Organism(name="Arabidopsis thaliana")
		athaliana.save()

	def _create_motifs(self):
		for comb in combinations:
			motif = Motif(motif=comb)
			motif.save()

	def _create_proteins(self):
		athaliana = Organism.objects.get(name="Arabidopsis thaliana")

		for rec in SeqIO.parse("Arabidopsis_filtered.fa", "fasta"):
			protein = Protein(
				organism = athaliana,
				proteinId = rec.id,
				description = rec.description,
				length = len(rec.seq),
				sequence = rec.seq
				)
			protein.save()

	def _create_results(self):
		for comb in combinations2:
			with open("fullresults/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
				jsonFile = json.load(f)

			for prot in jsonFile["proteins"]:
				result = Result(
					motif = Motif.objects.get(motif=comb),
					protein = Protein.objects.get(proteinId = prot["id"]),
					count = prot["count"],
					zscore = prot["zscore"],
					ratio = prot["ratio"],
					pvalue = prot["pvalue"]

					)
				result.save()

	def handle(self, *args, **options):
		self._create_results()



combinations2 = ['GW']

combinations = ['AA', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL',
'AM', 'AN', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AV', 'AW', 'AY', 'CC', 'CD', 'CE',
'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CP', 'CQ', 'CR', 'CS', 'CT',
'CV', 'CW', 'CY', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN',
'DP', 'DQ', 'DR', 'DS', 'DT', 'DV', 'DW', 'DY', 'EE', 'EF', 'EG', 'EH', 'EI',
'EK', 'EL', 'EM', 'EN', 'EP', 'EQ', 'ER', 'ES', 'ET', 'EV', 'EW', 'EY', 'FF',
'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP', 'FQ', 'FR', 'FS', 'FT', 'FV',
'FW', 'FY', 'GG', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN', 'GP', 'GQ', 'GR', 'GS',
'GT', 'GV', 'GW', 'GY', 'HH', 'HI', 'HK', 'HL', 'HM', 'HN', 'HP', 'HQ', 'HR',
'HS', 'HT', 'HV', 'HW', 'HY', 'II', 'IK', 'IL', 'IM', 'IN', 'IP', 'IQ', 'IR',
'IS', 'IT', 'IV', 'IW', 'IY', 'KK', 'KL', 'KM', 'KN', 'KP', 'KQ', 'KR', 'KS',
'KT', 'KV', 'KW', 'KY', 'LL', 'LM', 'LN', 'LP', 'LQ', 'LR', 'LS', 'LT', 'LV',
'LW', 'LY', 'MM', 'MN', 'MP', 'MQ', 'MR', 'MS', 'MT', 'MV', 'MW', 'MY', 'NN',
'NP', 'NQ', 'NR', 'NS', 'NT', 'NV', 'NW', 'NY', 'PP', 'PQ', 'PR', 'PS', 'PT',
'PV', 'PW', 'PY', 'QQ', 'QR', 'QS', 'QT', 'QV', 'QW', 'QY', 'RR', 'RS', 'RT',
'RV', 'RW', 'RY', 'SS', 'ST', 'SV', 'SW', 'SY', 'TT', 'TV', 'TW', 'TY', 'VV',
'VW', 'VY', 'WW', 'WY', 'YY']