from Bio import SeqIO
import subprocess
import os
import json
import csv

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

extractedOrthologs = {}
for key, value in fileLocations.items():
    with open(value[0]) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        extractedOrthologs[key] =[[row[1],row[2]] for row in reader]
     

combinations = ['VY','GI']

for comb in combinations:
    print("Calculating {}".format(comb))
    with open("best_zscores/{}-{}_bz.json".format(comb,comb[::-1]), 'r') as f:
        best = json.load(f)

    orthologs ={}
    for protein in best['proteins']:
        orthologs[protein["id"]] = []

    for key,value in extractedOrthologs.items():
        print("A.thaliana id \t zscore \t {} id \t length  count   zscore".format(key))
        for protein in best['proteins']:
            for row in value:
                if row[0] == protein["id"] and row[1]!="":
                    for rec in  SeqIO.parse(fileLocations[key][1], "fasta"):
                        if rec.id == row[1]:
                            orthologs[protein["id"]].append(rec)
                            count = rec.seq.count(comb) + rec.seq.count(comb[::-1])
                            print("{} \t {:.2f} \t {} \t {} \t {} \t {:.2f}".format(
                                protein["id"],
                                protein["zscore"],
                                rec.id,len(rec.seq),
                                count,
                                (count - best["averageOccurence"])/best["standardDeviation"])
                            )



    menu=''
    for protein in best['proteins']:
        menu = menu + '<a href="{0}_{1}-{2}.html" target="_blank">{2}</a>\n'.format(comb,comb[::-1],protein['id'])
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
                proteinDescription = rec.description
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

        linkToIndex = "<a href='{}_{}-index.html' target='_blank'>index</a>\n".format(comb,comb[::-1])

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
        <h3>{}</h3>
        <pre>
{}






{}

{}
        </pre>
    </body>
</html>
            '''.format(comb,comb[::-1],protein,proteinDescription,table,highlightedAlignment,linkToIndex)
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
