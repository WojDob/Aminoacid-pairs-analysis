import matplotlib.pyplot as plt
import json
import os
import numpy as np

for subdir, dirs, files in os.walk('./'):
    for file in files:
        with open(file, 'r') as f:
            try:
                motif = json.load(f)
                print(motif['motif'])
                ratios= [ x['ratio'] for x in motif['proteins'] if x['ratio']>1]

                plt.hist(ratios, bins = 50)
                plt.xlabel('Ratio')
                plt.ylabel('Number of proteins')
                locs, labels = plt.yticks()            
                plt.yticks(np.arange(0, 5001, step=250))  
                # locs, labels = plt.xticks()            
                # plt.xticks(np.arange(0, 25, step=1))            
                plt.title('{}'.format(motif['motif']))
                figure = plt.gcf() # get current figure
                figure.set_size_inches(15, 10)
                plt.savefig('./plots/{}.png'.format(motif['motif'][0]),dpi = 100)
                plt.clf()
                f.close()
            except:
                print('{} error'.format(motif['motif']))
