from Bio import SeqIO
import subprocess
import os
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

    orthologs ={}
    for protein in best['proteins']:
        orthologs[protein["id"]] = []


    print("A.thaliana id \t zscore \t\t B.rapa id \t length  count   zscore")
    for protein in best['proteins']:
        for row in brapa:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Brassica_rapa.Brapa_1.0.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        orthologs[protein["id"]].append(rec)
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))
    
    print("A.thaliana id \t zscore \t\t G.max id \t length  count   zscore")
    for protein in best['proteins']:
        for row in gmax:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Glycine_max.Glycine_max_v2.1.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        orthologs[protein["id"]].append(rec)
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))

    print("A.thaliana id \t zscore \t\t O.indica id \t length  count   zscore")
    for protein in best['proteins']:
        for row in indica:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Oryza_indica.ASM465v1.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        orthologs[protein["id"]].append(rec)
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))

    print("A.thaliana id \t zscore \t\t O.japonica id \t length  count   zscore")
    for protein in best['proteins']:
        for row in japonica:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Oryza_sativa.IRGSP-1.0.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        orthologs[protein["id"]].append(rec)
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))

    print("A.thaliana id \t zscore \t\t M.truncatula id \t length  count   zscore")
    for protein in best['proteins']:
        for row in mtrunc:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Medicago_truncatula.MedtrA17_4.0.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        orthologs[protein["id"]].append(rec)
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))

    print("A.thaliana id \t zscore \t\t P.patens id \t length  count   zscore")
    for protein in best['proteins']:
        for row in ppatens:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Physcomitrella_patens.Phypa_V3.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        orthologs[protein["id"]].append(rec)
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))

    print("A.thaliana id \t zscore \t\t V.vinifera id \t length  count   zscore")
    for protein in best['proteins']:
        for row in vvini:
            if row[0] == protein["id"] and row[1]!="":
                for rec in SeqIO.parse("data/Vitis_vinifera.12X.pep.all.fa", "fasta"):
                    if rec.id == row[1]:
                        orthologs[protein["id"]].append(rec)
                        count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                        print("{} \t {} \t {} \t {} \t {} \t {}".format(protein["id"],protein["zscore"],rec.id,len(rec.seq),count,(count - best["averageOccurence"])/best["standardDeviation"]))

    menu=''
    for protein in best['proteins']:
        menu = menu + '<a href="{0}_{1}-{2}.html">{2}</a>\n'.format(comb,comb[::-1],protein['id'])
        protein['table'] ='''length  {0}
count   {1}
ratio   {2:.2f}
zscore  {3:.2f}'''.format(protein['length'],protein['count'],protein['ratio'],protein['zscore'])

    for protein in orthologs:
        #write protein and its orthologs to file
        f = open("temp.txt","w")
        for rec in SeqIO.parse("data/Arabidopsis_filtered.fa", "fasta"):
            if rec.id == protein:
                SeqIO.write(rec, f, "fasta")
        for ortholog in orthologs[protein]:
            SeqIO.write(ortholog, f, "fasta")
        f.close()

        #run mafft and get the output
        mafftOutput = subprocess.check_output("mafft --auto --clustalout temp.txt", shell=True)

        #calculate where the motifs are in the output
        positions = list()
        for i in range(len(mafftOutput)-1):
            if mafftOutput[i]+mafftOutput[i+1] in [comb,comb[::-1]] :
                positions.append(i)
                positions.append(i+1)

        #replace motif letters with html
        mafftOutputList= list(mafftOutput)
        for i in range(len( mafftOutputList)):
            if i in positions:
                mafftOutputList[i] = '<span class="highlight">{}</span>'.format(mafftOutputList[i])

        highlightedAlignment = ''.join(mafftOutputList)

        for bestProtein in best['proteins']:
            if protein == bestProtein['id']:
                table  = bestProtein['table']

        htmlFile=open('alignments/{}_{}-{}.html'.format(comb,comb[::-1],protein),'w')
        content = '''
<!doctype html>
<html>
    <head>
        <title>{}-{}</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <link rel="stylesheet" href="styles.css">
    </head>
    <body>
        <h1>{}</h1>
        <pre>
{}






{}

{}
        </pre>
    </body>
</html>
            '''.format(comb,comb[::-1],protein,table,highlightedAlignment,menu)
        htmlFile.write(content)
        htmlFile.close()



    indexFile = open('alignments/{}_{}-index.html'.format(comb,comb[::-1]),"w")
    content = '''
<!doctype html>
<html>
    <head>
        <title>{0}-{1}</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <link rel="stylesheet" href="styles.css">
    </head>
    <body>
        <h1>{0}-{1}</h1>
        <pre>
{2}
        </pre>
    </body>
</html>
        '''.format(comb,comb[::-1],menu)
    indexFile.write(content)
    indexFile.close()
