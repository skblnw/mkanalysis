#!/usr/bin/env python

import os
import glob
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = [6, 4]
font = {'family': 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size': 16}
plt.rc('font', **font)
cmap = cm.get_cmap('Set1')


def customized_boxplot(parts):
    # cmap = ['g','#F6C344']
    # cmap = cm.get_cmap('Set1')
    for ii, pc in enumerate(parts['boxes']):
        pc.set_facecolor(cmap(ii))
        pc.set_edgecolor('none')
        pc.set_alpha(.5)

import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')
cmd = ap.parse_args()

prefix_list = ['free', 'bound']
trial_list = "trial[0-9]*"

data = []
for prefix in prefix_list:
    data_violin = np.empty(0, dtype=float)
    for trial in sorted(glob.glob(prefix+'/'+trial_list)):
        if Path(trial).is_dir():
            # data_profile = np.loadtxt(trial+"/fepout", comments=['#','@'])
            data_violin = np.append(data_violin, np.loadtxt(trial+"/fepout", comments=['#','@'])[-1,-1])
        else:
            print(trial+' does not have fepout!')
            quit()
    print('Sample Size: %d' % np.size(data_violin))
    print('Mean: %.2f' % np.mean(data_violin))
    print('s.d.: %.2f' % np.std(data_violin, ddof=1))
    data.append(data_violin)

print('')
print('ddG: %.2f' % ( np.mean(data[1])-np.mean(data[0]) ))
print('SE: %.2f' % ( np.sqrt(np.std(data[0], ddof=1)*np.std(data[0], ddof=1) / np.size(data[0]) + np.std(data[1], ddof=1)*np.std(data[1], ddof=1) / np.size(data[1]))))

fig, ax = plt.subplots()
parts = ax.boxplot(data,
    notch = False,
    showmeans = True,
    meanline = False,
    patch_artist=True,
    boxprops = dict(lw=2,facecolor='red'),
    whiskerprops = dict(lw=2, color='#7d758a'),
    capprops = dict(lw=2, color='#7d758a'),
    meanprops = dict(marker='D', markersize=8, markerfacecolor='#7d758a', markeredgecolor='none'),
    medianprops = dict(ls='-', lw=2, color='#7d758a'))

# ax.set_xlim([0,180])
# ax.set_ylim([274,286])
# ax.set_ylim([-274,-286])
# ax.set_xticks([])
# ax.set_yticks([])

# plt.xlabel(label)
# plt.ylabel(cmd.ylabel)
# plt.title(cmd.title)
# ax.set_yticks(np.arange(0,-400,-100))
ax.set_xticks([])
customized_boxplot(parts)

outfile = cmd.output
interactive = cmd.interactive

if interactive:
    plt.show()
else:
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight', pad_inches=.01)
