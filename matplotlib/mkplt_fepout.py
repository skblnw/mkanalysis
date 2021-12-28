#!/usr/bin/env python

import os
import glob
from pathlib import Path
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.font_manager
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


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, max(vals))
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, min(vals), q1)
    return lower_adjacent_value, upper_adjacent_value

def plot_profile(axes, data):
    axes.plot(data[:,1], data[:,-1], linewidth=2, marker="o")

def plot_window(axes, data):
    barplot = axes.bar(data[:,1], data[:,-2], width=0.01)
    for ii, pc in enumerate(barplot.patches):
        pc.set_edgecolor('none')
        pc.set_alpha(.3)

def plot_violin(axes, data):
    parts = axes.violinplot(data,
        showmeans=True, showmedians=False, showextrema=False)
    cmap = ['#DCDCDC','#D3772C']
    for ii, pc in enumerate(parts['bodies']):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('black')
        pc.set_alpha(.8)
    parts['cmeans'].set_color('#D43F3A')
    parts['cmeans'].set_lw(3)
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=0)
    whiskers = adjacent_values(data, quartile1, quartile3)
    whiskers_min, whiskers_max = whiskers[0], whiskers[1]
    axes.scatter(1, medians, marker='o', color='white', s=30, zorder=3)
    axes.vlines(1, quartile1, quartile3, color='k', linestyle='-', lw=5)
    axes.vlines(1, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    axes.set_xlim([.5, 1.5])
    axes.set_xticks([])
    axes.yaxis.set_label_position("right")
    axes.yaxis.tick_right()
    axes.title.set_text('%.2f +/- %.2f [%d]' % (data.mean(), (np.std(data, ddof=1)/np.sqrt(np.size(data))), np.size(data)))

import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', nargs='+', help='Path to fepout')
ap.add_argument('--window', help='Show dG of each window', action='store_true')
io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')
cmd = ap.parse_args()


prefix_trial = "trial[0-9]*"

fig = plt.figure(constrained_layout=True, figsize=(16, 8))
ax1, ax2 = fig.subplots(1,2, gridspec_kw={'width_ratios': [6, 1]})
data_violin = np.empty(0, dtype=float)
for path in cmd.input:
    for trial in sorted(glob.glob(path+'/'+prefix_trial)):
        if Path(trial).is_dir():
            data_profile = np.loadtxt(trial+"/fepout", comments=['#','@'])
            if cmd.window:
                plot_window(ax1, data_profile)
            else:
                plot_profile(ax1, data_profile)
            data_violin = np.append(data_violin, np.loadtxt(trial+"/fepout", comments=['#','@'])[-1,-1])
        else:
            print(trial+' does not have fepout!')
            quit()
# print(' Sample Size: %d' % np.size(data_violin))
# print(' Mean: %.2f' % data_violin.mean())
# print(' SEM: %.2f' % (np.std(data_violin, ddof=1)/np.sqrt(np.size(data_violin))))
plot_violin(ax2, data_violin)

# ax1.set_xlim([0,0.02])

if cmd.interactive:
    plt.show()
else:
    plt.tight_layout()
    plt.savefig(cmd.outfile)
