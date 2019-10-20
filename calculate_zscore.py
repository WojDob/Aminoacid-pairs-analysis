from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import sys
import json

frequency = {
    'A':0.082,
    'R':0.055,
    'N':0.04,
    'D':0.054,
    'C':0.013,
    'Q':0.039,
    'E':0.067,
    'G':0.07,
    'H':0.022,
    'I':0.059,
    'L':0.096,
    'K':0.058,
    'M':0.024,
    'F':0.038,
    'P':0.047,
    'S':0.066,
    'T':0.053,
    'W':0.01,
    'Y':0.029,
    'V':0.068
}


def pairFrequency(comb):
    a = list(comb)
    return frequency[a[0]]*frequency[a[1]]*2*100


comb = 'WG'

zScoreOutput = {
    'motif': [comb,comb[::-1]],
    'averageOccurence' : -1,
    'standardDeviation' : -1,
    'proteins':[]
}


countList = list()
for rec in SeqIO.parse("C:\\Users\\wojci\\Desktop\\github\\Aminoacid-pairs-analysis\\data\\Arabidopsis_filtered.fa", "fasta"):

    protein = {
        'id':rec.id,
        'length': len(rec.seq),
        'count':-1,
        'ratio':-1,
        'zscore':-1,
    }


    #calculate ratio of expected occurence and actual
    positions = list()
    for i in range(len(rec.seq)-1):
        if rec.seq[i]+rec.seq[i+1] in [comb,comb[::-1]] :
            positions.append(i)
            positions.append(i+1)
    percentage = len(set(positions))/len(rec.seq) *100

    protein['ratio'] = percentage / pairFrequency(comb)

    #count instances of motif
    recCombCount = rec.seq.count(comb) + rec.seq.count(comb[::-1])
    protein['count'] = recCombCount
    countList.append(recCombCount)

    zScoreOutput['proteins'].append(protein)

zScoreOutput["averageOccurence"] = np.mean(countList)
zScoreOutput["standardDeviation"] = np.std(countList)
#calculate zscore of each protein
for rec in zScoreOutput['proteins']:
    rec['zscore'] =(rec['count'] - zScoreOutput['averageOccurence'])/zScoreOutput["standardDeviation"]


#print(json.dumps(zScoreOutput,indent=4))
plt.hist([protein['zscore'] for protein in zScoreOutput['proteins']])
plt.show()
