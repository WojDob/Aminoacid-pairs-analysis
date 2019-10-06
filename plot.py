import matplotlib.pyplot as plt
import json
import os

for subdir, dirs, files in os.walk('./'):
    for file in files:
        with open(file, 'r') as f:
            try:
                motif = json.load(f)
                print(motif['motif'])
                ratios= [ x['ratio'] for x in motif['proteins']]

                plt.plot(ratios)
                plt.savefig('./plots/{}.png'.format(motif['motif'][0]),bbox_inches='tight')
                plt.clf()
                f.close()
            except:
                print("error {}".format(motif['motif']))
                pass
