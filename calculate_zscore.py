from Bio import SeqIO
import sys
import json
comb = 'WG'

zScoreOutput = {
    'motif': [comb,comb[::-1]],
    'averageOccurence' : -1,
    'standardDeviation' : -1,
    'proteins':[]
}



totalCombCount = 0
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    recCombCount = rec.seq.count(comb) + rec.seq.count(comb[::-1])
    protein = {
        'id':rec.id,
        'count':recCombCount,
        'zscore':-1,
    }
    totalCombCount += recCombCount
    zScoreOutput['proteins'].append(protein)

zScoreOutput["averageOccurence"] = totalCombCount / len(zScoreOutput["proteins"])

print(json.dumps(zScoreOutput,indent=4))
