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
