#!/usr/bin/env python

import os
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
    'axes.labelsize': 48,
    'axes.linewidth': 2,
    'axes.xmargin': 0,
    'lines.linewidth' : 3,
    'legend.fontsize': 24,
    'xtick.labelsize': 24,
    'xtick.major.size': 10,
    'xtick.major.width': 3,
    'ytick.labelsize': 40,
    'ytick.major.size': 10,
    'ytick.major.width': 3,
    'text.usetex': False,
    'figure.figsize': [16, 8],
    'figure.dpi': 100
   }
mpl.rcParams.update(params)


def customized_boxplot(parts):
    cmap = ['#505A5B','#D6FFF6','#6B4FCF','#4DCCBD','#2374AB']
    for ii, pc in enumerate(parts['boxes']):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('none')
        pc.set_alpha(.5)

def customized_bar(parts):
    cmap = ['#D6FFF6','#6B4FCF','#4DCCBD','#2374AB','#FF8484']
    cmap = ['#B0A990','#B0A990','#B0A990','#202030','#202030']
    for ii, pc in enumerate(parts.patches):
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

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, max(vals))
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, min(vals), q1)
    return lower_adjacent_value, upper_adjacent_value

def plot_profile(axes, data):
    axes.plot(data[:,1], data[:,-1], linewidth=2, marker="o")

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
    axes.title.set_text('n = %d' % np.size(data))

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


refpoints = {
    'L452R': [
            ['../../l452r/fepout_free']],
    'E484Q': [
            ['../../e484q/fepout_free']],
    'T478K': [
            ['../../t478k/fepout_free']],
    'L452R_E484Q': [
            ['../../l452r_e484q/fepout_free']],
    'L452R_T478K': [
            ['../../l452r_t478k/fepout_free']]
}

datapoints = {
    'L452R': [
            ['../../l452r/fepout_ace2']],
    'E484Q': [
            ['../../e484q/fepout_ace2']],
    'T478K': [
            ['../../t478k/fepout_ace2']],
    'L452R_E484Q': [
            ['../../l452r_e484q/fepout_ace2']],
    'L452R_T478K': [
            ['../../l452r_t478k/fepout_ace2']]
}

# /------------------/
# /     Figure 1     /
# /------------------/
ref = []
data = []
for points, title in zip([refpoints, datapoints], [ref, data]): 
    for label, values in points.items():
        print(label, values)
        fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw={'width_ratios': [4, 1]})
        data_violin = np.empty(0, dtype=float)
        for value in values:
            for trial in os.listdir(value[0]):
                data_profile = np.loadtxt(value[0]+'/'+str(trial), comments=['#','@'])
                plot_profile(ax1, data_profile)
                data_violin = np.append(data_violin, np.loadtxt(value[0]+'/'+str(trial), comments=['#','@'])[-1,-1])
        print(' Sample Size: %d' % np.size(data_violin))
        print(' Mean: %.2f' % data_violin.mean())
        print(' SEM: %.2f' % (np.std(data_violin, ddof=1)/np.sqrt(np.size(data_violin))))
        plot_violin(ax2, data_violin)
        title.append(data_violin)
        if not interactive:
            if title is ref:
                plt.tight_layout()
                plt.savefig(outfile+'/PLOT_fepout_'+label.lower()+'_free')
            else:
                plt.tight_layout()
                plt.savefig(outfile+'/PLOT_fepout_'+label.lower()+'_bound_ace2')

refmean = [ dd.mean() for dd in ref ]
refsem = [ np.std(dd, ddof=1) / np.sqrt(np.size(dd)) for dd in ref ]
datamean = [ dd.mean() for dd in data ]
datasem = [ np.std(dd, ddof=1) / np.sqrt(np.size(dd)) for dd in data ]

plotmean = [ dd - rr for dd,rr in zip(datamean,refmean) ]
plotsem = [ np.sqrt(dd*dd + rr*rr) for dd,rr in zip(datasem,refsem) ]
print('Mean: ')
for ii in plotmean: print(' %.2f' % ii) 
print('SEM: ')
for ii in plotsem: print(' %.2f' % ii) 

fig, ax = plt.subplots()
parts = ax.bar(np.arange(1,len(plotmean)+1),plotmean,
    yerr=plotsem,
    error_kw=dict(lw=2, capsize=12, capthick=2))
# plt.xlabel(cmd.xlabel)
# plt.ylabel(cmd.ylabel)
# plt.title(cmd.title)
# plt.ylim([-3,9.5])
# plt.yticks([0,5])
customized_bar(parts)
labels = ['L452R','E484Q','T478K','L452R_E484Q','L452R_T478K']
set_axis_style(ax, labels)
plt.xticks(rotation=60)

if interactive:
    plt.show()
else:
    plt.tight_layout()
    plt.savefig(outfile+'/PLOT_fepbox_ace2_ddG')
