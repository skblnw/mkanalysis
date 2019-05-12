#!/usr/bin/env python3

from pylab import *
import matplotlib.pyplot as plt


import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', type=str, help='Input', metavar='Input')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')


cmd = ap.parse_args()

params = {
   'axes.labelsize': 12,
   'font.size': 12,
   'legend.fontsize': 14,
   'xtick.labelsize': 20,
   'ytick.labelsize': 20,
   'text.usetex': False,
   'figure.figsize': [6, 6]
   }
rcParams.update(params)

data = np.loadtxt(cmd.input, comments=['#', '@'])
n_series = len(data[1,:])
n_elements = len(data[:,1])
print('[+] Read {0} series of data ({1} elements)'.format(n_series, n_elements))

axes(frameon=0)
f = plt.figure(1)
ax = plt.gca()
plt.grid(color="0.9", linestyle='-', linewidth=1)

ax.plot(data[:,0]/1000, data[:,1]*10, alpha=0.9, linewidth=2.5, linestyle='-', label='Complex')

#xlim(0, 300)
ylim(0, 6)
plt.xlabel('Time (ns)', fontsize=14)
plt.ylabel('RMSD (A)', fontsize=14)

legend = legend(loc=4);
frame = legend.get_frame()
frame.set_facecolor('0.9')
frame.set_edgecolor('0.9')


outfile=cmd.output
interactive=cmd.interactive
if outfile:
#    plt.tight_layout()
    plt.savefig(outfile)
        
if interactive:
    plt.show()
