#!/home/kevin/anaconda3/bin/python
import argparse
from argparse import RawDescriptionHelpFormatter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
params = {
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.linewidth': 2,
    'axes.xmargin': 0,
    'axes.ymargin': 0,
    'lines.linewidth' : 2,
    'legend.fontsize': 16,
    'xtick.labelsize': 16,
    'xtick.major.size': 8,
    'xtick.major.width': 3,
    'xtick.direction': 'in',
    'ytick.labelsize': 16,
    'ytick.major.size': 8,
    'ytick.major.width': 3,
    'ytick.direction': 'in',
    'errorbar.capsize': 20,
    'figure.figsize': [4, 4],
    'figure.dpi': 300
   }
mpl.rcParams.update(params)

# /----------------------------------/
# /         Argument Parser          /
# /----------------------------------/

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
ap.add_argument('data', type=str, help='Data for plotting')
ap.add_argument('--res', type=int, help='Residue Starting ID', nargs='?', default=1)
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')
io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')
parser = ap.parse_args()

# /----------------------------------/
# /            Read Data             /
# /----------------------------------/

resid_start = parser.res

f = open(parser.data, 'r')
firstline = f.readline()
size=int(firstline.split()[0])

data_raw = [line.strip().split() for line in f.readlines()]
data_raw[-1] = data_raw[-1][:-1]
data = [data for xx in data_raw for data in xx]

npa = np.array(data, dtype='float')
new = npa.reshape(size,size)

# /----------------------------------------------------/
# /                     Plotting Area                  /
# /----------------------------------------------------/

fig, ax = plt.subplots()
im = ax.imshow(new, 
                interpolation='gaussian', cmap=plt.cm.BuGn, 
                norm=colors.TwoSlopeNorm(vmin=0.3, vcenter=.5, vmax=1),
                origin='lower',
                extent=[resid_start,resid_start-1+size,resid_start,resid_start-1+size])

ax.set_xlim(430, 515)
ax.set_ylim(430, 515)
ax.set_xticks(range(445,515,15))
ax.set_yticks(range(445,515,15))

cbar = plt.colorbar(im,
                    ticks=[.4,.6,.8,1],
                    fraction=0.046, pad=0.04)
cbar.outline.set_linewidth(1)
cbar.ax.tick_params(length=2, width=1)

plt.tight_layout()
if parser.interactive:
    plt.show()
else:
    plt.savefig(parser.output)
