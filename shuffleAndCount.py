from Bio import SeqIO
import random
import utils

combinations = ['WG']

def countMotif(seq,comb):
	return seq.count(comb) + seq.count(comb[::-1])


allAminoacids = ""
proteinLengths = [0]

for rec in SeqIO.parse("data/test.fa", "fasta"):
	allAminoacids += rec.seq
	proteinLengths.append(proteinLengths[-1] + len(rec.seq))

allAminoacidsList = list(allAminoacids)

for i in range(1):
	print("Shuffling nr  " + str(i))
	random.shuffle(allAminoacidsList)

	fakeProteome = list()
	for i in range(len(proteinLengths)-1):
		fakeProteome.append(''.join(allAminoacidsList[proteinLengths[i] : proteinLengths[i+1]]))

	for comb in combinations:
		for sequence in fakeProteome:
			print(countMotif(sequence,comb))
	print(fakeProteome)
