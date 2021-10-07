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
    'font.size': 24,
    'axes.labelsize': 24,
    'axes.linewidth': 2,
    'axes.xmargin': 0,
    'lines.linewidth' : 3,
    'legend.fontsize': 24,
    'xtick.labelsize': 24,
    'xtick.major.size': 10,
    'xtick.major.width': 3,
    'ytick.labelsize': 36,
    'ytick.major.size': 10,
    'ytick.major.width': 3,
    'text.usetex': False,
    'figure.figsize': [16, 8],
    'figure.dpi': 100
   }
mpl.rcParams.update(params)


def customized_violinplot(parts):
    cmap = ['#353535','#3C6E71','#FFFFFF','#D9D9D9','#284B63']
    for ii, pc in enumerate(parts['bodies']):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('black')
        pc.set_alpha(.8)

    parts['cmeans'].set_color('#D43F3A')
    parts['cmeans'].set_lw(3)

def customized_boxplot(parts):
    cmap = ['#D6FFF6','#6B4FCF','#4DCCBD','#2374AB','#FF8484']
    for ii, pc in enumerate(parts['boxes']):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('none')
        pc.set_alpha(.5)

def add_whiskerbox(ax,data):
    quartiles = np.array([
        customized_percentile(column)
        for column in data])
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(data, quartiles[:,0], quartiles[:,2])])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(quartiles[:,1]) + 1)
    ax.scatter(inds, quartiles[:,1], marker='o', color='white', s=30, zorder=3)
    ax.vlines(inds, quartiles[:,0], quartiles[:,2], color='k', linestyle='-', lw=5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, max(vals))
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, min(vals), q1)
    return lower_adjacent_value, upper_adjacent_value

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


fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw={'width_ratios': [4, 1]})
ax1twin = ax1.twinx()
data_plot = np.empty(0, dtype=float)

# datapoints = list(range(1,4))+list(range(6,12))+list(range(13,18))+[19,20]
# datapoints = [1,2,3,11,12,13]
datapoints = range(1,31)

for ii in datapoints:

    filename = str(ii)

    tmp = np.loadtxt(filename, comments=['#','@'])
    data_plot = np.append(data_plot, tmp[-1,3])
    # ax1twin.bar(tmp[:,1], tmp[:,2])
    ax1.plot(tmp[:,1], tmp[:,3], linewidth=2, marker="o")

ax1twin.set_yticks([])

parts = ax2.violinplot(data_plot,
    showmeans=True, showmedians=False, showextrema=False)
cmap = ['#DCDCDC','#D3772C']
for ii, pc in enumerate(parts['bodies']):
    pc.set_facecolor(cmap[ii])
    pc.set_edgecolor('black')
    pc.set_alpha(.8)

parts['cmeans'].set_color('#D43F3A')
parts['cmeans'].set_lw(3)
quartile1, medians, quartile3 = np.percentile(data_plot, [25, 50, 75], axis=0)
whiskers = adjacent_values(data_plot, quartile1, quartile3)
whiskers_min, whiskers_max = whiskers[0], whiskers[1]
ax2.scatter(1, medians, marker='o', color='white', s=30, zorder=3)
ax2.vlines(1, quartile1, quartile3, color='k', linestyle='-', lw=5)
ax2.vlines(1, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

ax2.set_xlim([.5, 1.5])
ax2.set_xticks([])
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
ax2.title.set_text('n = %d' % np.size(data_plot))


datamean = data_plot.mean()
datasem = np.std(data_plot, ddof=1) / np.sqrt(np.size(data_plot))
print('Sample Size: %d' % np.size(data_plot))
print('Mean: %.2f' % datamean)
print('SEM: %.2f' % datasem)


# set style for the axes
# labels = ['1','2','3','1+3','4']
# set_axis_style(ax, labels)

outfile = cmd.output
interactive = cmd.interactive

if interactive:
    plt.show()
else:
    plt.tight_layout()
    plt.savefig(outfile)
