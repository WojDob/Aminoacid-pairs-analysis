from Bio import SeqIO
import random
import utils
import json
import decimal


def countMotif(seq,comb):
	if comb != comb[::-1]:
		return seq.count(comb) + seq.count(comb[::-1])
	else:
		return seq.count(comb)



numberOfCycles = 100
numberOfProteinsInProteome = 27628
proteome = "data/Arabidopsis_filtered.fa"
# proteome = "data/testShort.fa"
#combinations = ['WG','WW','AA','IW', 'IY','GG']
combinations = ['WG']

allAminoacids = ""
proteinLengths = [0]
motifCounts = {}

for comb in combinations:
	motifCounts[comb] = {}

#put all aminoacids from proteome into one string ang get the lengths of each protein 
for rec in SeqIO.parse(proteome, "fasta"):
	allAminoacids += rec.seq
	proteinLengths.append(proteinLengths[-1] + len(rec.seq))

allAminoacidsList = list(allAminoacids)

for i in range(numberOfCycles):

	print("Shuffling nr  " + str(i))
	random.shuffle(allAminoacidsList)

	fakeProteome = list()
	for i in range(len(proteinLengths)-1):
		fakeProteome.append(''.join(allAminoacidsList[proteinLengths[i] : proteinLengths[i+1]]))
	#print(fakeProteome)

	for comb in combinations:
		for sequence in fakeProteome:
			count = countMotif(sequence,comb)
			if count in motifCounts[comb]:
				motifCounts[comb][count] += 1
			else:
				motifCounts[comb][count] = 1

print(json.dumps(motifCounts, indent = 4))

best = {}
for comb in combinations:
	print(comb)
	with open("best_zscores/{}-{}_bz.json".format(comb,comb[::-1]), 'r') as f:
		best[comb] = json.load(f)
	for protein in best[comb]["proteins"]:
		sumOfProteinsWithEqualOrBiggerCount = 0
		for key in motifCounts[comb]:
			if key >= protein["count"]:
				sumOfProteinsWithEqualOrBiggerCount+=motifCounts[comb][key]
				
		print(str(sumOfProteinsWithEqualOrBiggerCount) + " " + str(protein["count"])+ " " + protein["id"])
		protein["pvalue"] = float(decimal.Decimal(sumOfProteinsWithEqualOrBiggerCount) / decimal.Decimal(numberOfProteinsInProteome * numberOfCycles))
print(json.dumps(best['WG'],indent = 4))