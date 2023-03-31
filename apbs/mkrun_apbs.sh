#!/bin/bash

[ $# -eq 0 ] && { echo "mkgmx> Usage: $0 [PDB]; Suggested: 0"; exit 1; }

PDB=$1

files=("$PDB")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

ln -s $PDB input.pdb

mkvmdfix input.pdb
pdb2pqr --ff AMBER input.pdb charged.pqr

cat > mg_auto.in<<EOF
read
    mol pqr charged.pqr
end
elec
    mg-auto
    dime 161 193 193
    cglen 98.8893 119.7177 123.182
    fglen 78.1702 90.4222 92.46
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 1.0
    sdie 78.54
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx elec-pot
end
quit
EOF

apbs mg_auto.in
