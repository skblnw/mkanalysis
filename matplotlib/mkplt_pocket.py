#!/usr/bin/env python

import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.font_manager

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
    'figure.figsize': [16, 16],
    'figure.dpi': 100
   }
matplotlib.rcParams.update(params)

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
# /      Figure 1: heatmap     /
# /---------------------------/

df = pd.read_csv("2", sep=" ")
del df['0']
cmap=["red","salmon","pink"]
sns.heatmap(df.T, cmap=cmap)


if interactive:
    plt.show()
