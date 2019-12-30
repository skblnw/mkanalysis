#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "DejaVu Sans"
import matplotlib.style
import matplotlib as mpl
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
mpl.style.use('publication')


import argparse
from argparse import RawDescriptionHelpFormatter

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
data = np.loadtxt(cmd.input, comments=['#','@'])

f = plt.figure(1)
ax = plt.gca()
ax.plot(data[:,0], data[:,1], linewidth=1)
plt.xlabel(cmd.xlabel)
plt.ylabel(cmd.ylabel)
plt.title(cmd.title)

outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.tight_layout()
    plt.savefig(outfile)
        
if interactive:
    plt.show()
