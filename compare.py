from Bio import SeqIO
import json
import csv

with open('data/hom_Brapa.txt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    brapa =[[row[1],row[2]] for row in reader]

with open('data/hom_Gmax.txt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    gmax =[[row[1],row[2]] for row in reader]

with open('data/hom_indica.txt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    indica =[[row[1],row[2]] for row in reader]

with open('data/hom_japonica.txt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    japonica =[[row[1],row[2]] for row in reader]

with open('data/hom_Mtranc.txt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    mtrunc =[[row[1],row[2]] for row in reader]

with open('data/hom_Ppatens.txt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    ppatens =[[row[1],row[2]] for row in reader]

with open('data/hom_Vvini.txt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    vvini =[[row[1],row[2]] for row in reader]


combinations = ['WG']

for comb in combinations:
    print(comb)
    with open("best_zscores/{}-{}_bz.json".format(comb,comb[::-1]), 'r') as f:
        best = json.load(f)
    print("A.thaliana id \t zscore \t\t B.rapa id \t length  count   zscore")
    for protein in best['proteins']:
        for row in brapa:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Brassica_rapa.Brapa_1.0.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))

    print("A.thaliana id \t zscore \t\t G.max id \t length  count   zscore")
    for protein in best['proteins']:
        for row in gmax:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Glycine_max.Glycine_max_v2.1.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))
    print("A.thaliana id \t zscore \t\t O.indica id \t length  count   zscore")
    for protein in best['proteins']:
        for row in indica:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Oryza_indica.ASM465v1.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))
    print("A.thaliana id \t zscore \t\t O.japonica id \t length  count   zscore")
    for protein in best['proteins']:
        for row in japonica:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Oryza_sativa.IRGSP-1.0.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))
    print("A.thaliana id \t zscore \t\t M.truncatula id \t length  count   zscore")
    for protein in best['proteins']:
        for row in mtrunc:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Medicago_truncatula.MedtrA17_4.0.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))
    print("A.thaliana id \t zscore \t\t P.patens id \t length  count   zscore")
    for protein in best['proteins']:
        for row in ppatens:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Physcomitrella_patens.Phypa_V3.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))
    print("A.thaliana id \t zscore \t\t V.vinifera id \t length  count   zscore")
    for protein in best['proteins']:
        for row in vvini:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Vitis_vinifera.12X.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))
