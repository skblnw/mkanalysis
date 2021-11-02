#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [4, 4]
font = {'family': 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size': 16}
plt.rc('font', **font)
import matplotlib.style
import matplotlib as mpl

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
for file in cmd.input:
    data = np.loadtxt(file, comments=['#','@'])
    # ax.plot(data[:,0], data[:,1], lw=.5, c='#3fc07d')
    ax.plot(data[:,0], data[:,1], lw=.5)

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
