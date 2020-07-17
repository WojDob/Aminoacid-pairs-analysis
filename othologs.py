from Bio import SeqIO
import subprocess
import os
import json
import csv


arabidopsisProtein = 'AT5G04290.1'
motif = 'GW'




fileLocations = {
    "B.rapa":["data/hom_Brapa.txt","data/Brassica_rapa.Brapa_1.0.pep.all.fa"],
    "G.max":["data/hom_Gmax.txt","data/Glycine_max.Glycine_max_v2.1.pep.all.fa"],
    "O.indica":["data/hom_indica.txt","data/Oryza_indica.ASM465v1.pep.all.fa"],
    "O.japonica":["data/hom_japonica.txt","data/Oryza_sativa.IRGSP-1.0.pep.all.fa"],
    "M.truncatula":["data/hom_Mtranc.txt","data/Medicago_truncatula.MedtrA17_4.0.pep.all.fa"],
    "P.patens":["data/hom_Ppatens.txt","data/Physcomitrella_patens.Phypa_V3.pep.all.fa"],
    "V.vinifera":["data/hom_Vvini.txt","data/Vitis_vinifera.12X.pep.all.fa"],
    "C.reingardtii":["data/hom_Creinhardtii.txt","data/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.pep.all.fa"],
    "Z.mays":["data/hom_Zmays.txt","data/Zea_mays.B73_RefGen_v4.pep.all.fa"],
}

def extractOrthologs(fileLocations):
    extractedOrthologs = {}
    for key, value in fileLocations.items():
        with open(value[0]) as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            extractedOrthologs[key] =[[row[1],row[2]] for row in reader]
    return extractedOrthologs

def findOrthologIdentifier(extractedOrthologs,identifier):
    orthologs = []
    for key in extractedOrthologs:
        for protein in extractedOrthologs[key]:
            if protein[0] == identifier and protein[1]!='' :
                print(key)
                print("{}: {}".format(protein[0],protein[1]))
                orthologs.append(protein[1])
    return orthologs
            


def findProteinAndCountPvalue(fileLocations,identifiers,comb):
    print("id   motif    count    pvalue\n")
    for identifier in identifiers:
        for key, value in fileLocations.items():
            for rec in SeqIO.parse(value[1], "fasta"):
                if rec.id == identifier:
                    #count instances of motif
                    if comb != comb[::-1]:
                        recCombCount = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                    else:
                        recCombCount = rec.seq.count(comb)

                    print("{} {} {} {}".format(
                        rec.id,
                        comb,
                        recCombCount,
                        calculatePvalue(recCombCount,comb)
                        ))
                    break

def calculatePvalue(count,comb):
    with open("fakeMotifCounts.json", 'r') as f:
        motifCounts = json.load(f)

    sumOfProteinsWithEqualOrBiggerCount = 0
    for key in motifCounts[comb]:
        if int(key) >= count:
            sumOfProteinsWithEqualOrBiggerCount+=motifCounts[comb][key]
    return sumOfProteinsWithEqualOrBiggerCount



findProteinAndCountPvalue(fileLocations,findOrthologIdentifier(extractOrthologs(fileLocations),arabidopsisProtein),motif)