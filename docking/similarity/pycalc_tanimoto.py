#!/bin/python3
# -Prepare combined.sdf: obabel *.mol2 -osdf > combined.sdf
# -Loop over mol2 in current directory 
# -Calculate tanimoto coefficient using OpenBabel
# -Print coeff>.6
# python3 script.py
# Kev Oct2018@OSU

from __future__ import print_function
import os
import subprocess

def read( ref, target_list ):
    return subprocess.check_output(['obabel', ref, target_list, '-ofpt'], stderr=subprocess.STDOUT).decode('utf-8')

def create_dict(raw, coeff):
    for line in raw.splitlines():
        if len(line.split()) is 6 and line.split()[0][1:] != title and float(line.split()[5]) > .6:
            coeff[line.split()[0][1:]] = line.split()[5]
    return

def print_dict(coeff):
    if coeff:
        print(title, end=': ')
        for key in sorted(coeff):
            print(key, end=' ')
        print("")

if __name__== "__main__":

    for filename in sorted(os.listdir(".")):
        if filename.endswith(".mol2"):
            raw=read(filename,'combined.sdf')
            title=raw.splitlines()[0][1:]
            
            coeff={}
            create_dict(raw, coeff)
            
            print_dict(coeff)
        else:
            continue
