#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 250
plt.rcParams['figure.figsize'] = [3,2]
font = {'family': 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size': 12}
plt.rc('font', **font)
# cmap = ['k','r','g','#F6C344']

import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', type=argparse.FileType('r'), nargs='+', help='Filename')
ap.add_argument('--scaling', type=float, help='Scaling factor', default=1)
ap.add_argument('--col', type=int, help='Select column')
ap.add_argument('--errcol', type=int, help='Select column')
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')
ap.add_argument('--colormap', type=str, help='Colormap', nargs='?', default='Set1')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')

cmd = ap.parse_args()
cmap = cm.get_cmap(cmd.colormap)

f = plt.figure(1)
ax = plt.gca()
for ii, file in enumerate(cmd.input):
    data = np.loadtxt(file, comments=['#','@'])
    if len(data.shape) == 1:
        ax.plot(data, lw=2, label=str(ii+1), c=cmap(ii))
    elif data.shape[1] < 3:
        ax.plot(data[:,-1], lw=2, label=str(ii+1), c=cmap(ii))
    elif cmd.col is not None:
        if cmd.errcol is not None: 
            print("Plotting column "+str(cmd.col))
            # print(np.min(data[:,cmd.col]))
            # plt.errorbar(data[:,0], data[:,cmd.col], yerr=data[:,cmd.errcol], lw=1, c=cmap(ii)) 
            ax.plot(data[:,cmd.col-1]*cmd.scaling, lw=1, c=cmap(ii)) 
            plt.fill_between(np.arange(len(data[:,cmd.col-1])), (data[:,cmd.col-1]-data[:,cmd.errcol-1])*cmd.scaling, (data[:,cmd.col-1]+data[:,cmd.errcol-1])*cmd.scaling, facecolor=cmap(ii), alpha=0.3)
        else:
            print("Plotting column "+str(cmd.col))
            ax.plot(data[:,cmd.col-1]*cmd.scaling, lw=1, label=str(ii), c=cmap(ii)) 
    else:
        print("File "+str(ii+1)+" :"+str(data.shape))
        for jj in range(1,6):
            print("Plotting column "+str(jj))
            ax.plot(data[:,jj], lw=2, label=str(jj), c=cmap(jj))

# ax.set_xlim([22,50])
# ax.set_xticks([30,40,50])
# ax.set_ylim([0,14])
# ax.set_yticks([0,1,2])

# plt.legend()
plt.xlabel(cmd.xlabel)
plt.ylabel(cmd.ylabel)
plt.title(cmd.title)

# plt.grid(axis='y')
# plt.axhline(y=41.12, ls='--', lw=1, color='k')

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    # plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight', pad_inches=.02)
        
if interactive:
    plt.show()
