#!/usr/bin/env python
from Bio import SeqIO
import json
import sys

combinations = ['AA', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL',
'AM', 'AN', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AV', 'AW', 'AY', 'CC', 'CD', 'CE',
'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CP', 'CQ', 'CR', 'CS', 'CT',
'CV', 'CW', 'CY', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN',
'DP', 'DQ', 'DR', 'DS', 'DT', 'DV', 'DW', 'DY', 'EE', 'EF', 'EG', 'EH', 'EI',
'EK', 'EL', 'EM', 'EN', 'EP', 'EQ', 'ER', 'ES', 'ET', 'EV', 'EW', 'EY', 'FF',
'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP', 'FQ', 'FR', 'FS', 'FT', 'FV',
'FW', 'FY', 'GG', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN', 'GP', 'GQ', 'GR', 'GS',
'GT', 'GV', 'GW', 'GY', 'HH', 'HI', 'HK', 'HL', 'HM', 'HN', 'HP', 'HQ', 'HR',
'HS', 'HT', 'HV', 'HW', 'HY', 'II', 'IK', 'IL', 'IM', 'IN', 'IP', 'IQ', 'IR',
'IS', 'IT', 'IV', 'IW', 'IY', 'KK', 'KL', 'KM', 'KN', 'KP', 'KQ', 'KR', 'KS',
'KT', 'KV', 'KW', 'KY', 'LL', 'LM', 'LN', 'LP', 'LQ', 'LR', 'LS', 'LT', 'LV',
'LW', 'LY', 'MM', 'MN', 'MP', 'MQ', 'MR', 'MS', 'MT', 'MV', 'MW', 'MY', 'NN',
'NP', 'NQ', 'NR', 'NS', 'NT', 'NV', 'NW', 'NY', 'PP', 'PQ', 'PR', 'PS', 'PT',
'PV', 'PW', 'PY', 'QQ', 'QR', 'QS', 'QT', 'QV', 'QW', 'QY', 'RR', 'RS', 'RT',
'RV', 'RW', 'RY', 'SS', 'ST', 'SV', 'SW', 'SY', 'TT', 'TV', 'TW', 'TY', 'VV',
'VW', 'VY', 'WW', 'WY', 'YY']

frequency = {
    'A':0.082,
    'R':0.055,
    'N':0.04,
    'D':0.054,
    'C':0.013,
    'Q':0.039,
    'E':0.067,
    'G':0.07,
    'H':0.022,
    'I':0.059,
    'L':0.096,
    'K':0.058,
    'M':0.024,
    'F':0.038,
    'P':0.047,
    'S':0.066,
    'T':0.053,
    'W':0.01,
    'Y':0.029,
    'V':0.068
}

def pairFrequency(comb):
    a = list(comb)
    return frequency[a[0]]*frequency[a[1]]*2*100


for comb in combinations:

    print("Calculating {}-{}".format(comb,comb[::-1]))
    f = open("output/{}-{}_ratio.json".format(comb,comb[::-1]),"w+")
    output = list()

    for rec in SeqIO.parse(sys.argv[1], "fasta"):
        positions = list()
        for i in range(len(rec.seq)-1):
            if rec.seq[i]+rec.seq[i+1] in [comb,comb[::-1]] :
                positions.append(i)
                positions.append(i+1)
        percentage = len(set(positions))/len(rec.seq) *100
        output.append((rec.id,len(rec.seq),percentage,pairFrequency(comb),percentage/pairFrequency(comb)))

    output.sort(key=lambda tup: tup[4], reverse = True)
    json.dump(output,f, indent =1)
    f.close()
