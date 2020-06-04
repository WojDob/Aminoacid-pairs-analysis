import json
import utils
import matplotlib.pyplot as plt

def countAminoacid(aminoacid, pvalueCounts):
    s = 0
    for key in pvalueCounts:
        if aminoacid in key:
            # print("{}  {} ".format(key,pvalueCounts[key]))
            s+= pvalueCounts[key]
    # print("sum of all {} motif proteins".format(aminoacid))
    print(aminoacid+aminoacid)
    print(pvalueCounts[aminoacid+aminoacid])
    return s


def plot_pvalue_of_proteins(boundary):

def plot_pvalue_of_proteins(boundary):
    pvalueCounts = {}
    allSiginificantProteins= []
    for comb in utils.combinations:

        with open("full_fixed_results/{}-{}.json".format(comb,comb[::-1]), 'r') as f:
            resultPage = json.load(f)

        for protein in resultPage["proteins"]:
            if protein["pvalue"] <= boundary:
                # if comb == 'GW':
                #    print(protein["id"])
                if comb in pvalueCounts:
                    pvalueCounts[comb]+=1
                else:
                    pvalueCounts[comb]=1
                protein["motif"] = comb
                allSiginificantProteins.append(protein)
    # print("sum of all proteins")
    # print(sum(pvalueCounts.values()))
    # print("number of motifs engaged")
    # print(len(pvalueCounts))

    # aminoacidDict = {}
    # for key in utils.frequency:
    #     aminoacidDict[key] = countAminoacid(key,pvalueCounts)
    # print(aminoacidDict)

    with open("allSignificantProteins.json", 'w') as f:
        json.dump(allSiginificantProteins, f ,indent = 4)



    # plt.bar(range(len(aminoacidDict)), list(aminoacidDict.values()), align='center')
    # plt.xticks(range(len(aminoacidDict)), list(aminoacidDict.keys()))
    # plt.xlabel('Aminokwas')
    # plt.ylabel('Suma znalezionych bialek z dipeptydow zawierajacych aminokwas')
    # plt.show() 

    # sump = 0
    # for key in pvalueCounts:
    #     if 'P' in key:
    #         print("{}  {} ".format(key,pvalueCounts[key]))
    #         sump+= pvalueCounts[key]
    # print("sum of all P motif proteins")
    # print(sump)
    # sumq = 0
    # for key in pvalueCounts:
    #     if 'Q' in key:
    #         print("{}  {} ".format(key,pvalueCounts[key]))
    #         sumq+= pvalueCounts[key]
    # print("sum of all Q motif proteins")
    # print(sumq)
    # sumg = 0
    # for key in pvalueCounts:
    #     if 'G' in key:
    #         print("{}  {} ".format(key,pvalueCounts[key]))
    #         sumg+= pvalueCounts[key]
    # print("sum of all G motif proteins")
    # print(sumg)

    # sume = 0
    # for key in pvalueCounts:
    #     if 'E' in key:
    #         print("{}  {} ".format(key,pvalueCounts[key]))
    #         sume+= pvalueCounts[key]
    # print("sum of all E motif proteins")
    # print(sume)

    # sums = 0
    # for key in pvalueCounts:
    #     if 'S' in key:
    #         print("{}  {} ".format(key,pvalueCounts[key]))
    #         sums+= pvalueCounts[key]
    # print("sum of all S motif proteins")
    # print(sums)

    # plt.bar(range(len(pvalueCounts)), list(pvalueCounts.values()), align='center')
    # plt.xticks(range(len(pvalueCounts)), list(pvalueCounts.keys()), fontsize=7)
    # plt.xlabel('Motyw')
    # plt.ylabel('Liczba bialek z p-value rownym lub mniejszym niz {}'.format(str(boundary)))
    # plt.show() 

    # for comb in utils.combinations:
    #     if comb not in pvalueCounts:
    #         print(comb)
    print(json.dumps(pvalueCounts, indent=4))


plot_pvalue_of_proteins(0.0001)