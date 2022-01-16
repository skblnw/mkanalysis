#!/bin/bash

for ii in revert
do
    echo "> Doing HOLE computation for $ii"
    rm -rf $ii; mkdir -p $ii; cd $ii
    cat >> hole.inp << EOF
coord ../$ii.pdb
radius ../amberuni.rad
cvect 0.0 0.0 1.0
sphpdb out.sph
endrad 30
EOF
    hole <hole.inp> out.txt
    egrep "mid-|sampled" out.txt | awk '{print $1" "$2}' > formatted.dat
    cd ..
done

echo "> Done. Results are named formatted.dat."
