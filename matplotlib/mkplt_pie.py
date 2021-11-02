#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = [4, 4]
font = {'family': 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size': 16}
plt.rc('font', **font)
import matplotlib.style
import matplotlib as mpl

import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

ap.add_argument('input', type=argparse.FileType('r'), nargs='+', help='Filename')
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')

io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')

cmd = ap.parse_args()
f = plt.figure(1)
ax = plt.gca()
data = []
for file in cmd.input:
    tmp = np.loadtxt(file, comments=['#','@'])
    # ax.plot(data[:,0], data[:,1], lw=.5, c='#3fc07d')
    data.append(tmp[:25])
data=np.concatenate(data)
prop = [np.mean(data, axis=0)[1:]]
prop = np.append(prop, 1-np.sum(np.mean(data,axis=0)[1:]))

print('Sample Size: %d' % data.shape[0])
print('H-bond: %.2f' % prop[0])
print('1 Water: %.2f' % prop[1])
print('2 Waters: %.2f' % prop[2])
print('No interaction: %.2f' % prop[3])

labels = ['H-bond', '1 Water', '2 Waters', 'No interaction']
explode = (0.1, 0.1, 0.1, 0)
colors = ['#8F91A2','#DCEDFF','#94B0DA','#505A5B']

wedges, texts = ax.pie(prop.tolist(), 
                                    colors=colors, 
                                    explode=explode, 
                                    wedgeprops=dict(width=0.5),
                                    # labels=labels, 
                                    # autopct='%1.1f%%', 
                                    shadow=False, 
                                    startangle=120)

bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    if prop[i] < 0.01:
        continue
    else:
        ax.annotate("{:.1%}".format(prop[i]), xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

# for text in texts:
#     text.set_color('grey')

plt.xlabel(cmd.xlabel)
plt.ylabel(cmd.ylabel)
plt.title(cmd.title)

outfile=cmd.output
interactive=cmd.interactive
if interactive:
    plt.show()
else:
    plt.tight_layout()
    plt.savefig(outfile)
