#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 150
plt.rcParams['figure.figsize'] = [4,4]
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
ap.add_argument('--xscaling', type=float, help='Scaling factor for x-axis', default=1)
ap.add_argument('--scaling', type=float, help='Scaling factor', default=1)
ap.add_argument('--yshift', type=float, help='Y-shift factor', default=0)
ap.add_argument('--col', type=int, nargs='+', help='Select column')
ap.add_argument('--errcol', type=int, help='Select column')
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')
ap.add_argument('--colormap', type=str, help='Colormap', nargs='?', default='Dark2')
ap.add_argument('--running-avg', action='store_true', help='Plot running averages')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')

cmd = ap.parse_args()
cmap = cm.get_cmap(cmd.colormap)

def moving_average(data, window_size=7):
    """
    Calculate the moving average of the given list or series.
    """
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

f = plt.figure(1)
ax = plt.gca()
for ii, file in enumerate(cmd.input):
    print(f"Plotting file: {file.name}")
    data = np.loadtxt(file, comments=['#','@'])
    if len(data.shape) == 1:
        if cmd.running_avg:
            ax.plot(data, lw=2, c=cmap(ii), alpha=.1)
            ax.plot(moving_average(data), label=f"File {ii}", c=cmap(ii))
        else:
            scaled_x = np.arange(len(data)) * cmd.xscaling
            ax.plot(scaled_x, data, lw=2, label=str(ii+1), c=cmap(ii))
    elif data.shape[1] < 3:
        if cmd.running_avg:
            ax.plot(data[:,1], color=cmap(ii), alpha=.2)
            ax.plot(moving_average(data[:,1]), lw=2, label=f"{ii+1}", color=cmap(ii))
            plt.legend(['HLA','HLA (groove)','b2M','peptide'])
        else:
            ax.plot(data[:,0]*cmd.xscaling, data[:,1], lw=2, label=str(ii+1), c=cmap(ii))
    elif cmd.col is not None:
        if cmd.errcol is not None: 
            col=cmd.col[0]
            print("Plotting column "+str(cmd.col))
            # print(np.min(data[:,cmd.col]))
            # plt.errorbar(data[:,0], data[:,cmd.col], yerr=data[:,cmd.errcol], lw=1, c=cmap(ii)) 
            ax.plot(data[:,0], data[:,col-1]*cmd.scaling+cmd.yshift, lw=1.5, c=cmap(ii)) 
            plt.fill_between(data[:,0], (data[:,col-1]-data[:,cmd.errcol-1])*cmd.scaling+cmd.yshift, (data[:,col-1]+data[:,cmd.errcol-1])*cmd.scaling+cmd.yshift, facecolor=cmap(ii), alpha=0.3)
            # plt.fill_between(data[:,0], np.arange(len(data[:,col-1])), (data[:,col-1]-data[:,cmd.errcol-1])*cmd.scaling, (data[:,col-1]+data[:,cmd.errcol-1])*cmd.scaling, facecolor=cmap(ii), alpha=0.3)
        else:
            for jj, col in enumerate(cmd.col):
                print("Plotting column "+str(col))
                ax.plot(data[:,0]*cmd.xscaling, data[:,col-1]*cmd.scaling, lw=1, label=str(col), c=cmap(jj)) 
    else:
        print("File "+str(ii+1)+" :"+str(data.shape))
        for jj in range(1,data.shape[1]):
            print("Plotting column "+str(jj))
            ax.plot(data[:,0]*cmd.xscaling, data[:,jj], lw=2, label=str(jj), c=cmap(jj))

# ax.set_xlim([33,48])
# ax.set_xticks([])
# ax.set_ylim([0,8])
# ax.set_yticks([0,1,2,3])

# ax.set_ylim([0,21])
# ax.set_yticks([0,10,20])

ax.set_xlim(left=0)
# ax.set_ylim([0,10])
# ax.set_yticks([0,-2,-4,-6,-8,-10])

def offset_formatter(x, pos):
    return f'{int(x + 0)}'
ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(offset_formatter))

# plt.legend(['C-ter','N-ter'])
# legend = plt.legend(loc='best')
# custom_labels = ['HLA','HLA (groove)','b2M','peptide']
# for text, label in zip(legend.get_texts(), custom_labels):
#     text.set_text(label)
# plt.legend(['State Entropy'])

plt.xlabel(cmd.xlabel)
plt.ylabel(cmd.ylabel)
plt.title(cmd.title)

# plt.grid(axis='y')
# plt.axhline(y=0, ls='-', lw=1, color='k')

f.tight_layout()

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.savefig(outfile, bbox_inches='tight', pad_inches=.02)
        
if interactive:
    plt.show()
