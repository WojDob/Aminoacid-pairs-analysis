import json
import decimal



def countMotif(seq,comb):
	if comb != comb[::-1]:
		return seq.count(comb) + seq.count(comb[::-1])
	else:
		return seq.count(comb)

def pvalue(sum,numOfProteins,numOfCycles):
	return float(decimal.Decimal(sum) / decimal.Decimal(27628 * 1000))

dic = {
        "0": 16095711, 
        "1": 7394000, 
        "2": 2749207, 
        "3": 925661, 
        "4": 303551, 
        "5": 100621, 
        "6": 34911, 
        "7": 13468, 
        "8": 5642, 
        "9": 2460, 
        "10": 1315, 
        "11": 703, 
        "12": 361, 
        "13": 184, 
        "14": 113, 
        "15": 44, 
        "16": 22, 
        "17": 15, 
        "18": 6, 
        "19": 2, 
        "20": 3
    }


with open("json_files/GW-WG.json", 'r') as f:
	WG = json.load(f)

for protein in WG["proteins"]:
	sumOfProteinsWithEqualOrBiggerCount = 0
	for key in dic:
		if int(key) >= protein["count"]:
			sumOfProteinsWithEqualOrBiggerCount+=dic[key]
	protein["pvalue"] = pvalue(sumOfProteinsWithEqualOrBiggerCount,27628,1000)

	
with open('GW-WG.json','w') as f:
        json.dump(WG, f ,indent = 4)