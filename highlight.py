from Bio import SeqIO
import utils

identifier = 'AT1G67870.1'
comb = 'HQ'

for rec in SeqIO.parse("data/Arabidopsis_filtered.fa", "fasta"):
    if rec.id == identifier:
        positions = list()
        for i in range(len(rec.seq)-1):
            if rec.seq[i]+rec.seq[i+1] in [comb,comb[::-1]] :
                positions.append(i)
                positions.append(i+1)

        #replace motif letters with html
        recSeqList= list(rec.seq)
        print(len(rec.seq))
        for i in range(len( recSeqList)):
            if i in positions:
                recSeqList[i] = '<span class="highlight">{}</span>'.format(recSeqList[i])
        i = 60
        while i < len(recSeqList):
            recSeqList.insert(i, '<br>')
            i += (60+1)
        highlightedMotifs = ''.join(recSeqList)

        content = '''
<!doctype html>
<html>
    <head>
       
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <link rel="stylesheet" href="styles.css">
    </head>
    <body>

        <pre>
{}

        </pre>
    </body>
</html>
            '''.format(highlightedMotifs)



with open('{}_{}-{}.html'.format(identifier,comb,comb[::-1]),'w') as htmlFile:
    htmlFile.write(content)
    # htmlFile.close()
