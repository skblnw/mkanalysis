#/usr/bin/python2
# -Prepare a namelist, the namelist could contain either ZINC# or # only
# -Download mol2 from ZINC website and generate png using obabel
# python2 script.py
# Kev Oct2018@OSU

import re
import urllib
import subprocess

f = open('namelist', 'r')
for line in f:
    index = line.strip()
    print index
    url = "http://zinc.docking.org/substance/"+index
    f = urllib.urlopen(url)

    for line in f.readlines():
        found = re.findall(r'Download MOL2 File', line)
        if found:
            line = line.split('"')[1]
            urllib.urlretrieve(line, index+".mol2")
            cmd = "obabel "+index+".mol2 -O "+index+".png -d"
            subprocess.Popen(cmd, shell=True)
            break
