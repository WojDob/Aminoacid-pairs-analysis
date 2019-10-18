from Bio import SeqIO
import numpy as np
import sys
import json
comb = 'WG'

zScoreOutput = {
    'motif': [comb,comb[::-1]],
    'averageOccurence' : -1,
    'averageFrequency' : -1,
    'standardDeviation' : -1,
    'proteins':[]
}



totalCombCount = 0
countList = list()
for rec in SeqIO.parse("C:\\Users\\wojci\\Desktop\\github\\Aminoacid-pairs-analysis\\data\\test.fa", "fasta"):
    recCombCount = rec.seq.count(comb) + rec.seq.count(comb[::-1])
    protein = {
        'id':rec.id,
        'count':recCombCount,
        'ratio':-1,
        'zscore':-1,
    }
    countList.append(recCombCount)
    totalCombCount += recCombCount
    zScoreOutput['proteins'].append(protein)

zScoreOutput["averageOccurence"] = np.mean(countList)
zScoreOutput["standardDeviation"] = np.std(countList)


print(json.dumps(zScoreOutput,indent=4))
print(countList)
