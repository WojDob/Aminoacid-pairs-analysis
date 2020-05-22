import json
import utils
import matplotlib.pyplot as plt


def plot_pvalue_of_proteins(boundary):
    pvalueCounts = {}
    for comb in utils.combinations:

        with open("full_fixed_results/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
            resultPage = json.load(f)

        for protein in resultPage["proteins"]:
            if protein["pvalue"] <= boundary:
                if comb in pvalueCounts:
                    pvalueCounts[comb]+=1
                else:
                    pvalueCounts[comb]=1
    print(sum(pvalueCounts.values()))
    print(len(pvalueCounts))

    print(pvalueCounts['PP'])
    sump = 0
    for key in pvalueCounts:
        if 'P' in key:
            print("{}  {} ".format(key,pvalueCounts[key]))
            sump+= pvalueCounts[key]
    print(sump)

    plt.bar(range(len(pvalueCounts)), list(pvalueCounts.values()), align='center')
    plt.xticks(range(len(pvalueCounts)), list(pvalueCounts.keys()), fontsize=7)
    plt.xlabel('motif')
    plt.ylabel('Liczba bialek z p-value rownym lub mniejszym niz {}'.format(str(boundary)))
    plt.show() 


plot_pvalue_of_proteins(0.00001)