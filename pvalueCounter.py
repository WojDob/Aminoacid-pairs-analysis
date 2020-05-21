import json
import utils


combinations = ['AC']
pvalueCounts ={
    "0.01":0,
    "0.001":0,
    "0.0001":0,
    "0.00001":0,
    "0.000001":0,
    "0.0" : 0,

}
for comb in utils.combinations:
    with open("fullresults/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
        resultPage = json.load(f)
    for protein in resultPage["proteins"]:
        for key in pvalueCounts:
            if protein["pvalue"] <= float(key):
                pvalueCounts[key]+=1

print(pvalueCounts)