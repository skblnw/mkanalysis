#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl
import seaborn as sns
import scipy.stats as st
# plt.rcParams['axes.linewidth'] = 1
# plt.rcParams['figure.dpi'] = 100
# plt.rcParams['figure.figsize'] = [16, 4]
# font = {'family': 'sans-serif',
#         'sans-serif': 'Arial',
#         'style': 'normal',
#         'weight': 'normal',
#         'size': 16}
# plt.rc('font', **font)
params = {
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 48,
    'axes.labelsize': 48,
    'axes.linewidth': 2,
    'axes.xmargin': 0,
    'lines.linewidth' : 3,
    'legend.fontsize': 26,
    'xtick.labelsize': 40,
    'xtick.major.size': 10,
    'xtick.major.width': 3,
    'ytick.labelsize': 40,
    'ytick.major.size': 10,
    'ytick.major.width': 3,
    'text.usetex': False,
    'figure.figsize': [8, 8],
    'figure.dpi': 100
   }
mpl.rcParams.update(params)

def customized_histogram(parts):
    cmap = ['#7e817f','#35CAAD','#EA5A15']
    for ii, pc in enumerate(parts[2]):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('none')
        pc.set_alpha(.5)

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
colors = ['#7e817f','#35CAAD','#EA5A15']
linecolors = ['#463F3A', '#0C64F3', '#C44900']
for ii, file in enumerate(cmd.input):
    data = np.loadtxt(file, comments=['#','@'])

    q25, q75 = np.percentile(abs(data[:,1]), [0.25, 0.75])
    bin_width = 2 * (q75 - q25) * len(data) ** (-1/3)
    bins = round((abs(data[:,1]).max() - abs(data[:,1]).min()) / bin_width)
    print("Freedman–Diaconis number of bins:", bins)

    ax.hist(abs(data[:,1]), color=colors[ii], alpha = 0.95, density=True, bins=bins)
    print('Sample Size: %d' % np.size(data[:,1]))
    print('Mean: %.2f' % abs(data[:,1]).mean())
    print('SEM: %.2f' % np.std(abs(data[:,1]), ddof=1))

    mn, mx = plt.xlim()
    plt.xlim(mn, mx)
    kde_xs = np.linspace(mn, mx, 300)
    kde = st.gaussian_kde(abs(data[:,1]))
    ax.plot(kde_xs, kde.pdf(kde_xs), lw=3, color=linecolors[ii], label="PDF")

    # sns.displot(data[:,1], bins=bins, kde=True)

    ax.plot([851, 851], [0, .01], 'k--', lw=2)

plt.xlabel(cmd.xlabel)
plt.ylabel(cmd.ylabel)
plt.title(cmd.title)
ax.set_ylim([0, .01])
ax.set_yticks([0, .01])
ax.set_xlim([690, 1030])
ax.set_xticks([800, 1000])

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.tight_layout()
    plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=.01)
        
if interactive:
    plt.show()
