#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [4, 4]
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

ap.add_argument('--xfile', type=str, nargs='+', help='Filename')
ap.add_argument('--yfile', type=str, nargs='+', help='Filename')
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
for ii, (xfile, yfile) in enumerate(zip(cmd.xfile[0].split(), cmd.yfile[0].split())):
    print(ii, xfile, yfile)
    data1 = np.loadtxt(xfile, comments=['#','@'])
    data2 = np.loadtxt(yfile, comments=['#','@'])
    print("Shape of data1: "+str(data1.shape))
    print("Shape of data2: "+str(data2.shape))
    ax.scatter(data1[20:], data2[20:],color=cmap(ii))

# # rank1
# ax.scatter(205, 14.89, marker='*', color='k')
# # rank2
# ax.scatter(232, 52.45, marker='*', color='k')
# # 6dfx
# ax.scatter(247, 24.18, marker='*', color='r')
# # songyi
# ax.scatter(223, 10.97, marker='*', color='r')

ax.set_xlim([0,250])
ax.set_ylim([0,90])
# ax.set_xticks([])
# ax.set_yticks([])

plt.xlabel(cmd.xlabel)
plt.ylabel(cmd.ylabel)
plt.title(cmd.title)

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight', pad_inches=.01)
        
if interactive:
    plt.show()
