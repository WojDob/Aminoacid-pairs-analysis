from Bio import SeqIO
import json
import utils
top10 =[]
for i in range(15):
	topInIteration = {"count": 35, "length": 1024, "ratio": 0.0, "zscore": 38.63544976398659, "id": "AT3G28550.1"}
	for comb in utils.combinations:
	    print("Calculating {}".format(comb))

	    # with open("best_zscores/{}-{}_bz.json".format(comb,comb[::-1]), 'r') as f:
	    #     records = json.load(f)
	    
	    with open("full_fixed_results/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
	        records = json.load(f)	    

	    for protein in records['proteins']:
	    	protein['motif'] = comb
	    	if protein['count'] > topInIteration['count'] and protein not in top10:

	    		topInIteration = protein

	top10.append(topInIteration)
print(json.dumps(top10,indent = 4))