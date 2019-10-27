from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import sys
# import json
import utils


def pairFrequency(comb):
    a = list(comb)
    return utils.frequency[a[0]]*utils.frequency[a[1]]*2*100

def getDataList(key):
    return [protein[key] for protein in zScoreOutput['proteins']]



for comb in utils.combinations:

    zScoreOutput = {
        'motif': [comb,comb[::-1]],
        'averageOccurence' : -1,
        'standardDeviation' : -1,
        'proteins':[]
    }
    print("Calculating {}".format(zScoreOutput['motif']))
    countList = list()
    #for every protein record in data file
    for rec in SeqIO.parse(sys.argv[1], "fasta"):

        protein = {
            'id':rec.id,
            'length': len(rec.seq),
            'count':-1,
            'ratio':-1,
            'zscore':-1,
        }

        #calculate ratio of expected occurence and actual
        # positions = list()
        # for i in range(len(rec.seq)-1):
        #     if rec.seq[i]+rec.seq[i+1] in [comb,comb[::-1]] :
        #         positions.append(i)
        #         positions.append(i+1)
        # percentage = len(set(positions))/len(rec.seq) *100
        # protein['ratio'] = percentage / pairFrequency(comb)

        #count instances of motif
        recCombCount = rec.seq.count(comb) + rec.seq.count(comb[::-1])
        protein['count'] = recCombCount
        countList.append(recCombCount)

        #add the protein to the list
        zScoreOutput['proteins'].append(protein)

    zScoreOutput["averageOccurence"] = np.mean(countList)
    zScoreOutput["standardDeviation"] = np.std(countList)

    #calculate zscore of each protein
    for rec in zScoreOutput['proteins']:
        rec['zscore'] =(rec['count'] - zScoreOutput['averageOccurence'])/zScoreOutput["standardDeviation"]


    zScoresList = getDataList('zscore')
    # ratioList = getDataList('ratio')
    lengthList = getDataList('length')


    plt.scatter(zScoresList, lengthList,alpha=0.1)
    plt.xlabel('Z-score')
    plt.ylabel('Length of sequence')
    plt.title(zScoreOutput['motif'])
    figure = plt.gcf()
    figure.set_size_inches(15, 10)
    plt.savefig('./zscore_plots/zscore-length_scatterplots/{}-{}_z-l_scttr.png'.format(zScoreOutput['motif'][0],zScoreOutput['motif'][1]),dpi = 100)
    plt.clf()


    #
    # #make zscore histogram
    # plt.hist(zScoresList,bins = 50)
    # plt.title(zScoreOutput['motif'])
    # figure = plt.gcf()
    # figure.set_size_inches(15, 10)
    # plt.savefig('./zscore_plots/zscore_histograms/{}-{}_zs_hist.png'.format(zScoreOutput['motif'][0],zScoreOutput['motif'][1]),dpi = 100)
    # plt.clf()
    #
    # #make ratio to length scatterplot
    # plt.scatter(ratioList,lengthList,alpha=0.1)
    # plt.xlabel('Ratio')
    # plt.ylabel('Length of sequence')
    # plt.title(zScoreOutput['motif'])
    # figure = plt.gcf()
    # figure.set_size_inches(15, 10)
    # plt.savefig('./zscore_plots/ratio-length_scatterplots/{}-{}_r-l_scttr.png'.format(zScoreOutput['motif'][0],zScoreOutput['motif'][1]),dpi = 100)
    # plt.clf()
    #
    # #make ratio to zscore scatterplot
    # plt.scatter(ratioList,zScoresList,alpha=0.1)
    # plt.xlabel('Ratio')
    # plt.ylabel('Z-score')
    # plt.title(zScoreOutput['motif'])
    # figure = plt.gcf()
    # figure.set_size_inches(15, 10)
    # plt.savefig('./zscore_plots/ratio-zscore_scatterplots/{}-{}_r-zs_scttr.png'.format(zScoreOutput['motif'][0],zScoreOutput['motif'][1]),dpi = 100)
    # plt.clf()
