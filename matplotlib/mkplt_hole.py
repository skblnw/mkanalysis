#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 200
plt.rcParams['figure.figsize'] = [2, 4]
font = {'family': 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size': 16}
plt.rc('font', **font)
cmap = cm.get_cmap('Set2')

import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', type=argparse.FileType('r'), nargs='+', help='Filename')
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')

cmd = ap.parse_args()
f = plt.figure(1)
ax = plt.gca()
for ii, file in enumerate(cmd.input):
    data = np.loadtxt(file, comments=['#','@'])
    ax.plot(data[:,1], -data[:,0], lw=1.5, c=cmap(ii))
    ax.plot(-data[:,1], -data[:,0], lw=1.5, c=cmap(ii))

plt.axvline(x=-1, color='k', lw=.5, linestyle='--')
plt.axvline(x=1, color='k', lw=.5, linestyle='--')

ax.set_xlim([-6,6])
ax.set_ylim([-34,6])
# ax.set_ylim([-220,-180])
# ax.set_xticks([])
ax.set_yticks([])

plt.xlabel(cmd.xlabel)
plt.ylabel(cmd.ylabel)
plt.title(cmd.title)

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.tight_layout()
    plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=.01)
        
if interactive:
    plt.show()
