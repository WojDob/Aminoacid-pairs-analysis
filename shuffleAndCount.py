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

def pvalue(sum,numOfProteins,numOfCycles):
	return float(decimal.Decimal(sum) / decimal.Decimal(numOfProteins * numOfCycles))



numberOfCycles = 1000
numberOfProteinsInProteome = 27628
proteome = "data/Arabidopsis_filtered.fa"
combinations = utils.combinations


allAminoacids = ""
proteinLengths = [0]
motifCounts = {}

for comb in combinations:
	motifCounts[comb] = {}

print("Getting all aminoacids into a list")
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

	for comb in combinations:
		for sequence in fakeProteome:
			count = countMotif(sequence,comb)
			if count in motifCounts[comb]:
				motifCounts[comb][count] += 1
			else:
				motifCounts[comb][count] = 1

print(json.dumps(motifCounts, indent = 4))
with open('fakeMotifCounts.json','w') as f:
	json.dump(motifCounts,f ,indent =4)

print("Calculating results")
best = {}
for comb in combinations:
	print(comb)

	with open("json_files/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
		best[comb] = json.load(f)

	for protein in best[comb]["proteins"]:
		sumOfProteinsWithEqualOrBiggerCount = 0
		for key in motifCounts[comb]:
			if key >= protein["count"]:
				sumOfProteinsWithEqualOrBiggerCount+=motifCounts[comb][key]
		protein["pvalue"] = pvalue(sumOfProteinsWithEqualOrBiggerCount,numberOfProteinsInProteome,numberOfCycles)
		
	with open('fullresults/{}-{}.json'.format(comb,comb[::-1]),'w') as f:
	        json.dump(best[comb], f ,indent = 4)