#!/usr/bin/env python

import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
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
# matplotlib.style.use('tableau-colorblind10')
params = {
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 12,
    'axes.labelsize': 24,
    'axes.linewidth': 2,
    'axes.xmargin': 0,
    'lines.linewidth' : 1,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'xtick.major.size': 2,
    'xtick.major.width': 2,
    'ytick.labelsize': 12,
    'ytick.major.size': 2,
    'ytick.major.width': 2,
    'text.usetex': False,
    'figure.figsize': [8, 8],
    'figure.dpi': 100
   }
matplotlib.rcParams.update(params)

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, max(vals))
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, min(vals), q1)
    return lower_adjacent_value, upper_adjacent_value

def plot_profile(axes, data):
    axes.plot(data[:,1], data[:,-1], linewidth=1, marker=".")

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

ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')
cmd = ap.parse_args()
outfile = cmd.output
interactive = cmd.interactive

# /---------------------------/
# /      Figure 1: fepout     /
# /---------------------------/

prefix_posi = "p[0-9]*"
prefix_trial = "t[0-9]*"
states = ['free','bound']
df_mean = pd.DataFrame(columns=states)
df_sem = pd.DataFrame(columns=states)
dict_mean = {}
dict_sem = {}

fig = plt.figure(constrained_layout=True, figsize=(16, 12))
subfigs = fig.subfigures(len(sorted(glob.glob(prefix_posi))), 3, width_ratios=[1,8,8])

for ii, posi in enumerate(sorted(glob.glob(prefix_posi))): 
    datamean = []
    datasem = []
    subfigs[ii][0].text(0.5, 0.5, posi, fontsize=24, va="center", ha="center")
    for jj, state in enumerate(states): 
        # print('[Position '+posi+'] ['+state+' State]')

        ax1, ax2 = subfigs[ii][jj+1].subplots(1,2, gridspec_kw={'width_ratios': [6, 1]})

        data_violin = np.empty(0, dtype=float)
        for trial in sorted(glob.glob(posi+"/"+state+"/"+prefix_trial)):
            if Path(trial).is_dir():
                data_profile = np.loadtxt(trial+"/fepout", comments=['#','@'])
                plot_profile(ax1, data_profile)
                data_violin = np.append(data_violin, np.loadtxt(trial+"/fepout", comments=['#','@'])[-1,-1])
            else:
                print(trial+' does not have fepout!')
                quit()
        # print(' Sample Size: %d' % np.size(data_violin))
        # print(' Mean: %.2f' % data_violin.mean())
        # print(' SEM: %.2f' % (np.std(data_violin, ddof=1)/np.sqrt(np.size(data_violin))))
        plot_violin(ax2, data_violin)
        dict_mean[state] = data_violin.mean()
        dict_sem[state] = np.std(data_violin, ddof=1)/np.sqrt(np.size(data_violin))
    df_mean.loc[posi] = pd.Series(dict_mean)
    df_sem.loc[posi] = pd.Series(dict_sem)

if not interactive:
    plt.savefig(outfile+'_alascan_fepout')

# /------------------------/
# /      Figure 2: ddG     /
# /------------------------/
params = {
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 12,
    'axes.labelsize': 24,
    'axes.linewidth': 1,
    'axes.xmargin': 0,
    'lines.linewidth': 1,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'xtick.major.size': 4,
    'xtick.major.width': 1,
    'ytick.labelsize': 12,
    'ytick.major.size': 4,
    'ytick.major.width': 1,
    'figure.dpi': 300
   }
matplotlib.rcParams.update(params)

means = df_mean['bound'] - df_mean['free']
errors = np.sqrt(df_sem['bound']*df_sem['bound'] + df_sem['free']*df_sem['free'])

data=pd.read_csv("bb1", sep=" ", index_col=0, names=sorted(glob.glob(prefix_posi)))

fig = plt.figure(constrained_layout=True, figsize=(6, 3))
ax = means.plot.bar(color='grey', yerr=errors, capsize=3, rot=0)
ax.axhline(y=0, color='k')
ax2 = ax.twinx()
data.mean().plot(color='green', alpha=.6, yerr=data.std(), capsize=2, capthick=3)
data.mean().plot.area(color='green', alpha=.1)
ax2.set_ylim([1.2,.3])
# ax2.set_yticks(np.arange(0,1,0.2))
# ax2.set_yticklabels(np.arange(100,0,-20).tolist())


if interactive:
    plt.show()
else:
    plt.savefig(outfile+'_alascan_ddG')
