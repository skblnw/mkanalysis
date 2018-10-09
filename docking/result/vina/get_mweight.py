#/usr/bin/python2
# -Prepare a namelist, the namelist could contain either ZINC# or # only
# -Download mol2 from ZINC website and generate png using obabel
# python2 script.py
# Kev Oct2018@OSU

import re
import urllib
import subprocess
import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
ap.add_argument('input', type=str, help='Contains ZINC ID line by line', metavar='namelist')
cmd = ap.parse_args()

print "Name          Mol.W (g/mol)   xlogP"
f = open(cmd.input, 'r')
for line in f:
    index = line.strip()
    url = "http://zinc.docking.org/substance/"+index
    f = urllib.urlopen(url)

    lines = f.readlines()
    for ii in range(0, len(lines)):
        line = lines[ii]
        found = re.findall(r'td class=\"weight\" title=\"Molecular weight \(g/mol\)\"', line)
        if found:
            mweight=lines[ii+1].split()[0]
        found = re.findall(r'td class=\"xlogp\" title=\"Partition coefficient\"', line)
        if found:
            xlogp=lines[ii+1].split()[0]
        try:
            mweight
            xlogp
        except NameError:
            continue
        else:
            if float(xlogp) > 5:
                print index+"  "+mweight+"  "+xlogp+" XXX rejected due to high xlogp"
            elif float(xlogp) < 3:
                print index+"  "+mweight+"  "+xlogp+" OOO recommended due to low xlogp"
            else:
                print index+"  "+mweight+"  "+xlogp
            del mweight
            del xlogp
            break
