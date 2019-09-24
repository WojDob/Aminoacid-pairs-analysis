#!/usr/bin/env python
from Bio import SeqIO
import sys

combinations = ['AC', 'AD',
'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL', 'AM', 'AN', 'AP', 'AQ',
'AR', 'AS', 'AT', 'AV', 'AW', 'AY', 'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CK',
'CL', 'CM', 'CN', 'CP', 'CQ', 'CR', 'CS', 'CT', 'CV', 'CW', 'CY', 'DE', 'DF',
'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN', 'DP', 'DQ', 'DR', 'DS', 'DT', 'DV',
'DW', 'DY', 'EF', 'EG', 'EH', 'EI', 'EK', 'EL', 'EM', 'EN', 'EP', 'EQ', 'ER',
'ES', 'ET', 'EV', 'EW', 'EY', 'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP',
'FQ', 'FR', 'FS', 'FT', 'FV', 'FW', 'FY', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN',
'GP', 'GQ', 'GR', 'GS', 'GT', 'GV', 'GW', 'GY', 'HI', 'HK', 'HL', 'HM', 'HN',
'HP', 'HQ', 'HR', 'HS', 'HT', 'HV', 'HW', 'HY', 'IK', 'IL', 'IM', 'IN', 'IP',
'IQ', 'IR', 'IS', 'IT', 'IV', 'IW', 'IY', 'KL', 'KM', 'KN', 'KP', 'KQ', 'KR',
'KS', 'KT', 'KV', 'KW', 'KY', 'LM', 'LN', 'LP', 'LQ', 'LR', 'LS', 'LT', 'LV',
'LW', 'LY', 'MN', 'MP', 'MQ', 'MR', 'MS', 'MT', 'MV', 'MW', 'MY', 'NP', 'NQ',
'NR', 'NS', 'NT', 'NV', 'NW', 'NY', 'PQ', 'PR', 'PS', 'PT', 'PV', 'PW', 'PY',
'QR', 'QS', 'QT', 'QV', 'QW', 'QY', 'RS', 'RT', 'RV', 'RW', 'RY', 'ST', 'SV',
'SW', 'SY', 'TV', 'TW', 'TY', 'VW', 'VY', 'WY']

comb = combinations[1]

for rec in SeqIO.parse(sys.argv[1], "fasta"):
    print(rec.seq)
    positions = list()
    for i in range(len(rec.seq)-1):
        if rec.seq[i]+rec.seq[i+1] in [comb,comb[::-1]] :
            positions.append(i)
            positions.append(i+1)
    print(len(set(positions))/len(rec.seq))
