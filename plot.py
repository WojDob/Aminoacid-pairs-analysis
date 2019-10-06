import matplotlib.pyplot as plt
import json
import os

#kind of broken but not really, move it to output folder for it to even start

for subdir, dirs, files in os.walk('./'):
    for file in files:
        with open(file, 'r') as f:
            # try:
            motif = json.load(f)
            print(motif['motif'])
            ratios= [ x['ratio'] for x in motif['proteins']]

            plt.plot(ratios)
            plt.xlabel('Proteins')
            plt.ylabel('Ratio')
            plt.title('{} plot'.format(motif['motif']))
            plt.savefig('./plots/{}.png'.format(motif['motif'][0]),bbox_inches='tight')
            plt.clf()
            f.close()
            # except:
            #     print("error {}".format(motif['motif']))
