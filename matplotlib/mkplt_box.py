#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.font_manager
"""
'seaborn-darkgrid', 'Solarize_Light2', 'seaborn-notebook', 
'classic', 'seaborn-ticks', 'grayscale', 
'bmh', 'seaborn-talk', 'dark_background', 
'ggplot', 'fivethirtyeight', '_classic_test', 
'seaborn-colorblind', 'seaborn-deep', 'seaborn-whitegrid', 
'publication', 'seaborn', 'seaborn-poster', 
'seaborn-bright', 'seaborn-muted', 'seaborn-paper', 
'seaborn-white', 'fast', 'seaborn-pastel', 
'seaborn-dark', 'tableau-colorblind10', 'seaborn-dark-palette'
"""
"""
'Solarize_Light2', '_classic_test_patch', 'bmh', 'classic', 
'dark_background', 'fast', 'fivethirtyeight', 'ggplot', 'grayscale', 
'seaborn', 'seaborn-bright', 'seaborn-colorblind', 'seaborn-dark', 
'seaborn-dark-palette', 'seaborn-darkgrid', 'seaborn-deep', 'seaborn-muted', 
'seaborn-notebook', 'seaborn-paper', 'seaborn-pastel', 'seaborn-poster', 
'seaborn-talk', 'seaborn-ticks', 'seaborn-white', 'seaborn-whitegrid', 
'tableau-colorblind10'
"""
# mpl.style.use('tableau-colorblind10')
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
    'figure.figsize': [12, 12],
    'figure.dpi': 100
   }
mpl.rcParams.update(params)


def customized_boxplot(parts):
    cmap = ['#7e817f','#35CAAD','#EA5A15']
    for ii, pc in enumerate(parts['boxes']):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('none')
        pc.set_alpha(.5)

def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='in')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_tick_params(direction='in')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)


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



datapoints = {
    'WT': [
            ['rbdace2_zn',2]],
    'L452R+E484Q': [
            ['rbdace2_zn_l452r_e484q',2]],
    'L452R+T478K': [
            ['rbdace2_zn_l452r_t478k',2]]
}

data = []
for label, values in datapoints.items():
    print(label, values)
    tmp = np.empty(0, dtype=float)
    for value in values:
        name = value[0]
        ntrial = value[1]
        for trial in range(1, ntrial+1):
            tmp = np.append(tmp, np.loadtxt(name+'/t'+str(trial)+'.energy.xvg', comments=['#','@'])[:,1])
    print('Sample Size: %d' % np.size(tmp))
    print('Mean: %.2f' % tmp.mean())
    print('SEM: %.2f' % np.std(tmp, ddof=1))
    data.append(tmp)

fig, ax = plt.subplots()
# parts = ax.violinplot(data,
#     showmeans=True, showmedians=False, showextrema=False)
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
# plt.xlabel(cmd.xlabel)
# plt.ylabel(cmd.ylabel)
# plt.title(cmd.title)
# ax.set_yticks(np.arange(0,-400,-100))
ax.set_xticks([])
customized_boxplot(parts)

# set style for the axes
# labels = datapoints.keys()
# set_axis_style(ax, labels)

outfile = cmd.output
interactive = cmd.interactive

if interactive:
    plt.show()
else:
    plt.tight_layout()
    plt.savefig(outfile)
