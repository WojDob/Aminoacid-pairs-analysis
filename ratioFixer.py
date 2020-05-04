import json
import utils

combinations = ['GW']

for comb in combinations:
    print(comb)
    with open("fullresults/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
        fullresults = json.load(f)
    with open("old/ratio_data/{}-{}_ratio.json".format(comb,comb[::-1]), 'r') as f:
        ratioFile = json.load(f)

    for protein in fullresults["proteins"]:
    	for proteinWithRatio in ratioFile["proteins"]:
    		if protein["id"] == proteinWithRatio["id"]:
    			protein["ratio"] = proteinWithRatio["ratio"]
    			break


    with open('full_fixed_results/{}-{}.json'.format(comb,comb[::-1]),'w') as f:
        json.dump(fullresults, f ,indent = 4)