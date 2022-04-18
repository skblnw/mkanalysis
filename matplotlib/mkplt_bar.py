#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 200
plt.rcParams['figure.figsize'] = [2, 2]
font = {'family': 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size': 32}
plt.rc('font', **font)
cmap = cm.get_cmap('Set2')
# cmap = ['k','r','g','#F6C344']

import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', type=argparse.FileType('r'), nargs='+', help='Filename')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')

cmd = ap.parse_args()
# #354fca
fig = plt.figure(constrained_layout=True, figsize=(9, 3))
for ii, file in enumerate(cmd.input):
    data=pd.read_csv(file, sep=" ", index_col=0, names=['Region 1','Region 2A','Region 2B','Region 3'])
    ax = data[20:].mean().plot.bar(color=cmap(ii), yerr=data[20:].sem(), stacked=False, position=-ii+2, capsize=3, rot=0, width=.1)
    # ax = data[20:].mean().plot.bar(color=cmap(np.arange(len(data))), yerr=data[20:].sem(), stacked=False, position=ii+.5, capsize=3, rot=0, width=.5)
    # ax = data[20:].mean().plot.bar(color='grey', yerr=data[20:].sem(), stacked=False, position=ii+.5, capsize=3, rot=0, width=.5)
# plt.legend(['Region 1','Region 2A','Region 2B','Region 3'])
ax.set_xlim([-.5,3.5])
ax.set_ylim([0,10])
# ax.set_xticks([])
ax.set_yticks([0,2,4,6,8])

plt.grid(axis='y', lw=.3)

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.savefig(outfile, bbox_inches='tight', pad_inches=.01)
        
if interactive:
    plt.show()
