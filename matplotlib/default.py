#!/usr/bin/env python3

from pylab import *
import matplotlib.pyplot as plt


import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', type=str, help='Records of payment', metavar='Records of payment')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')


cmd = ap.parse_args()

data = np.loadtxt(cmd.input)

f = plt.figure(1)
ax = plt.gca()
ax.plot(data[:][0], data[:][1], '^')




outfile=cmd.output
interactive=cmd.interactive
if outfile:
    plt.tight_layout()
    plt.savefig(outfile)
        
if interactive:
    plt.show()