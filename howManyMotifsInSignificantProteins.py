import json
import utils
import matplotlib.pyplot as plt

with open("allSignificantProteins.json", 'r') as f:
    significantProteins = json.load(f)
idsOfAllProteins = [protein["id"] for protein in significantProteins]
motifCounts = dict((i, idsOfAllProteins.count(i)) for i in idsOfAllProteins)
print(json.dumps(motifCounts,indent=4))

countOfCounts = {}
for count in motifCounts:
    if motifCounts[count] in countOfCounts:
        countOfCounts[motifCounts[count]]+=1
    else:
        countOfCounts[motifCounts[count]]=1

print(json.dumps(countOfCounts,indent = 4))
plt.bar(range(len(countOfCounts)), list(countOfCounts.values()), align='center')
plt.xticks(range(len(countOfCounts)), list(countOfCounts.keys()))
plt.xlabel('Liczba nadreprezentowanych motywow')
plt.ylabel('Ilosc bialek')
plt.show() 
