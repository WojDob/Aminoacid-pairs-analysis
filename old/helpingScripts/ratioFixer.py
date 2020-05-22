import json
import utils
import decimal


def pvalue(sum,numOfProteins,numOfCycles):
    return float(decimal.Decimal(sum) / decimal.Decimal(numOfProteins * numOfCycles))

numberOfCycles = 1000
numberOfProteinsInProteome = 27628

with open("fakeMotifCounts.json", 'r') as f:
        motifCounts = json.load(f)


for comb in utils.combinations:
    print(comb)
    with open("fullresults/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
        fullresults = json.load(f)
    with open("old/ratio_data/{}-{}_ratio.json".format(comb,comb[::-1]), 'r') as f:
        ratioFile = json.load(f)

    for protein in fullresults["proteins"]:
        for proteinWithRatio in ratioFile["proteins"]:
            if protein["id"] == proteinWithRatio["id"]:
                protein["ratio"] = proteinWithRatio["ratio"]

                if comb == comb[::-1]:
                    protein["count"] = protein["count"]/2

                    sumOfProteinsWithEqualOrBiggerCount = 0
                    for key in motifCounts[comb]:
                        if int(key) >= protein["count"]:
                            sumOfProteinsWithEqualOrBiggerCount+=motifCounts[comb][key]
                    protein["pvalue"] = pvalue(sumOfProteinsWithEqualOrBiggerCount,numberOfProteinsInProteome,numberOfCycles)  

                break





    with open('full_fixed_results/{}-{}.json'.format(comb,comb[::-1]),'w') as f:
        json.dump(fullresults, f ,indent = 4)