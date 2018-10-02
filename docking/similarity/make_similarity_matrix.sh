#!/bin/bash
# -Create similarity matrices using OpenBabel
# Kev Oct2018@OSU

for ii in $(ls *.mol2);
do
    obabel $ii combined.sdf -ofpt | grep "^>" | tail -n +2 | cut -d= -f2 | awk 'BEGIN { ORS = " " } { print } END { printf( "\n" ); }' >> sim_fp2
    obabel $ii combined.sdf -ofpt -xfFP3 | grep "^>" | tail -n +2 | cut -d= -f2 | awk 'BEGIN { ORS = " " } { print } END { printf( "\n" ); }' >> sim_fp3
    obabel $ii combined.sdf -ofpt -xfFP4 | grep "^>" | tail -n +2 | cut -d= -f2 | awk 'BEGIN { ORS = " " } { print } END { printf( "\n" ); }' >> sim_fp4
    obabel $ii combined.sdf -ofpt -xfMACCS | grep "^>" | tail -n +2 | cut -d= -f2 | awk 'BEGIN { ORS = " " } { print } END { printf( "\n" ); }' >> sim_maccs
done
