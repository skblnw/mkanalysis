#!/home/kevin/anaconda3/bin/python
import argparse
from argparse import RawDescriptionHelpFormatter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.font_manager
from matplotlib import colors
params = {
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.linewidth': 1,
    'axes.xmargin': 0,
    'axes.ymargin': 0,
    'lines.linewidth' : 2,
    'legend.fontsize': 16,
    'xtick.labelsize': 12,
    'xtick.major.size': 4,
    'xtick.major.width': 1,
    'xtick.direction': 'in',
    'ytick.labelsize': 12,
    'ytick.major.size': 4,
    'ytick.major.width': 1,
    'ytick.direction': 'in',
    'errorbar.capsize': 20,
    'figure.figsize': [2, 4],
    'figure.dpi': 300
   }
mpl.rcParams.update(params)


def customized_percentile(data):
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=0)
    return quartile1, medians, quartile3

def customized_violinplot(parts):
    cmap = ['#353535','#3C6E71','#FFFFFF','#D9D9D9','#284B63']
    for ii, pc in enumerate(parts['bodies']):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('black')
        pc.set_alpha(.8)

    parts['cmeans'].set_color('#D43F3A')
    parts['cmeans'].set_lw(3)

def customized_boxplot(parts):
    cmap = ['#7e817f','#D6FFF6','#6B4FCF','#4DCCBD','#FF8484']
    for ii, pc in enumerate(parts['boxes']):
        pc.set_facecolor(cmap[ii])
        pc.set_edgecolor('none')
        pc.set_alpha(.5)

def customized_bar(parts):
    cmap = ['#D6FFF6','#6B4FCF','#4DCCBD','#FF8484']
    for ii, pc in enumerate(parts.patches):
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
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='in')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_tick_params(direction='in')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)


# /----------------------------------/
# /         Argument Parser          /
# /----------------------------------/

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', type=str, help='Filename')
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')
cmd = ap.parse_args()

# /----------------------------------/
# /            Read Data             /
# /----------------------------------/

reference = ['/data/kevin/sarscov2/7bwj/md/rbdace2/analysis/countHeavy/1.dat',
             '/data/kevin/sarscov2/7bwj/md/rbdace2/analysis/countHeavy/2.dat']
reference = ['result-wt-heavy.dat']
farray = [np.loadtxt(f, comments=['#','@']) for f in reference]
dataref = np.concatenate(farray)
# dataref = np.loadtxt(reference[1], comments=['#','@'])

data = np.loadtxt(cmd.input, comments=['#','@'])
dataplot = np.true_divide(np.subtract(data, dataref.mean(axis=0)), np.where(dataref.std(axis=0)<1, 1, dataref.std(axis=0)))
np.nan_to_num(dataplot, nan=0)

# /----------------------------------------------------/
# /                     Plotting Area                  /
# /----------------------------------------------------/

# create discrete colormap
# cmap = colors.ListedColormap(['white', '#677E98',  '#338FCC', '#0016FF'])
cmap = colors.ListedColormap(['#59A65E', 'white',  'white', '#A659A1'])
bounds = [-10, -1, 1, 10]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots()
im = ax.imshow(dataplot.transpose(), interpolation='gaussian', 
                                    cmap=cmap, norm=norm, 
                                    origin='lower', extent=[1,data.shape[0], 97,96+data.shape[1]],
                                    aspect=2
                                    )

cbar = plt.colorbar(im,
                    ticks=[],
                    fraction=0.046, pad=0.04)
cbar.outline.set_linewidth(1)
cbar.ax.tick_params(length=0, width=1)

# draw gridlines
# ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
# ax.set_yticks(np.arange(333, 526, 100))
# ax.set_yticks(np.arange(0, 40, 1));

# ax.set_ylim(430, 515)
# ax.set_yticks(range(445,515,15))
ax.set_xticks([])

plt.tight_layout()
interactive = cmd.interactive
if interactive:
    plt.show()
else:
    plt.savefig(cmd.output)
