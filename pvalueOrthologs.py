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
                #print(key)
                #print("{}: {}".format(protein[0],protein[1]))
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

def findSequence(identifier):
    if identifier.startswith("Bra"):
        for rec in SeqIO.parse("data/Brassica_rapa.Brapa_1.0.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec
    elif identifier.startswith("PNW"):
        for rec in SeqIO.parse("data/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec
    elif identifier.startswith("KR") or identifier.startswith("RC"):
        for rec in SeqIO.parse("data/Glycine_max.Glycine_max_v2.1.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec
    elif identifier.startswith("KEH") or identifier.startswith("AE"):
        for rec in SeqIO.parse("data/Medicago_truncatula.MedtrA17_4.0.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec
    elif identifier.startswith("BGI"):
        for rec in SeqIO.parse("data/Oryza_indica.ASM465v1.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec
    elif identifier.startswith("Os"):
        for rec in SeqIO.parse("data/Oryza_sativa.IRGSP-1.0.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec    
    elif identifier.startswith("Pp") or identifier.startswith("PAC"):
        for rec in SeqIO.parse("data/Physcomitrella_patens.Phypa_V3.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec
    elif identifier.startswith("VIT"):
        for rec in SeqIO.parse("data/Vitis_vinifera.12X.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec
    elif identifier.startswith("Zm"):
        for rec in SeqIO.parse("data/Zea_mays.B73_RefGen_v4.pep.all.fa", "fasta"):
            if rec.id == identifier:
                return rec    
    else:
        print("BUG {}".format(identifier))
        return False


#findProteinAndCountPvalue(fileLocations,findOrthologIdentifier(extractOrthologs(fileLocations),arabidopsisProtein),motif)


extractedOrthologs = extractOrthologs(fileLocations)

with open("allSignificantProteins.json", 'r') as f:
    allSignificantProteins = json.load(f)


protIds = [prot["id"] for prot in allSignificantProteins]

print(len(set(protIds)))

noOrthologs = 0
for prot in protIds:
    if len(findOrthologIdentifier(extractedOrthologs,prot)) == 0:
        noOrthologs+=1
print(noOrthologs)




# sanityCounter=0
# sumOfOrthologs = 0
# sumOfFound = 0
# for protein in allSignificantProteins:
#     sanityCounter+=1
#     print(sanityCounter)
#     #print("**{} {}**".format(protein["id"],protein["motif"]))
#     orthologs = findOrthologIdentifier(extractedOrthologs,protein["id"])
#     sumOfFound += len(orthologs)
#     #print(orthologs)
#     pacLock = False
#     for ortholog in orthologs:
#         if ortholog.startswith('PAC'):
#             pacLock = True;
#     if pacLock:
#         for ortholog in orthologs:
#             #print(ortholog)
#             rec = findSequence(ortholog)
#             if rec:
#                 comb = protein["motif"]
#                 if comb != comb[::-1]:
#                     recCombCount = rec.seq.count(comb) + rec.seq.count(comb[::-1])
#                 else:
#                     recCombCount = rec.seq.count(comb)
#                 if recCombCount > 4:
#                     if calculatePvalue(recCombCount,comb) <= 0.0001:
#                         sumOfOrthologs+=1
#                         break
# print(sumOfFound)
# print (sumOfOrthologs)
