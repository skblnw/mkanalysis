#!/usr/bin/env python

import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import colors
import matplotlib.style
import matplotlib.font_manager
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [4, 4]
font = {'family': 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size': 16}
plt.rc('font', **font)
cmap = cm.get_cmap('Reds_r')

import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

# ap.add_argument('input', type=argparse.FileType('r'), nargs='+', help='Filename')
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')

cmd = ap.parse_args()

prefix_posi = "smd1"
df = pd.DataFrame()

f = plt.figure(1)
ax = plt.gca()

for ii, file in enumerate(sorted(glob.glob(prefix_posi+"/*.dat"), key=lambda name: int(name[5:-4]))): 
    frame = np.loadtxt(file, comments=['#','@'])
    # ax.plot(frame[:,0], frame[:,1], lw=.2, color=cmap(ii), alpha=.1)
    # ax.scatter([ii]*len(frame[:,0]), frame[:,0], color=cmap(int(frame[:,1])), alpha=.1, s=10)
    sc = ax.scatter([ii]*len(frame[:,0]), frame[:,0], c=frame[:,1], cmap=cmap, vmin=0, vmax=3, alpha=1, s=99, marker='s')
    # for jj, posi in enumerate(frame[:,0]):
    #     ax.scatter(ii, posi, color=cmap(int(frame[jj,1])), alpha=.1, s=10)
    # ax.plot(data[:,0], data[:,1], lw=.5)
    # ax.plot(data, lw=.5)

ax.set_yticks([])
# ax.set_xlim([0,30])
# cbar = plt.colorbar(sc,
#                     ticks=[],
#                     fraction=0.046, pad=0.04)
# cbar.outline.set_linewidth(1)
# cbar.ax.tick_params(length=0, width=1)
# plt.colorbar(cm.ScalarMappable(norm=norm, cmap=frame[:,1]), ax=ax)
plt.colorbar(sc)

# plt.xlabel(cmd.xlabel)
# plt.ylabel(cmd.ylabel)
# plt.title(cmd.title)

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.tight_layout()
    plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=.01)
        
if interactive:
    plt.show()
