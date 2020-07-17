import json
import utils
import matplotlib.pyplot as plt

wantedCombinations = []

for comb in utils.combinations:
	if 'P' in comb or 'G' in comb or 'H' in comb or 'E' in comb:
		wantedCombinations.append(comb)


with open("allSignificantProteins.json", 'r') as f:
    significantProteins = json.load(f)
idsOfAllProteins = [protein["id"] for protein in significantProteins]

idsOfProteinsFromCombinations = []
for comb in wantedCombinations:
	for protein in significantProteins:
		if protein['motif'] == comb:
			idsOfProteinsFromCombinations.append(protein["id"])

print(wantedCombinations)
print(len(set(idsOfProteinsFromCombinations)))


setOfAllIds = set(idsOfAllProteins)

print(len(setOfAllIds))

# motifCounts = dict((i, idsOfAllProteins.count(i)) for i in idsOfAllProteins)
# print(json.dumps(motifCounts,indent=4))

# countOfCounts = {}
# for count in motifCounts:
#     if motifCounts[count] in countOfCounts:
#         countOfCounts[motifCounts[count]]+=1
#     else:
#         countOfCounts[motifCounts[count]]=1

# print(json.dumps(countOfCounts,indent = 4))
# plt.bar(range(len(countOfCounts)), list(countOfCounts.values()), align='center')
# plt.xticks(range(len(countOfCounts)), list(countOfCounts.keys()))
# plt.xlabel('Liczba nadreprezentowanych motywow')
# plt.ylabel('Liczba bialek')
# plt.show() 
