#!/bin/bash
# Kev @ 2019
# Bash script to calculate RMSD profile

[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [PREFIX] [REF]"; exit 1; }

PDB="$1"
TRJ="$2"
PREFIX="$3"
REF="$4"

SELREF="noh segname PROA"
REGIONS=("noh segname PROA" "noh segname PROB" "noh segname PROC" "index 58 59 96 97 99 101 102 104 106 140 142 693 910 912 915 976 977 978 979 980 1021 1023 1026 1029 1032 1035 1039 1041 1043 1045 1051 1083 1087 1088 1089 1091 1093 1096 1098 1099 1101 1102 1136 1138 1140 1142 1146 1147 1193 1195 1207 1209 1212 1213 1214 1216 1251 1264 1266 1270 1319 1320 1324 1514 1517 1549 1551 1552 1556 1784 1786 1815 1816 1818 1820 1821 1823 1825 1925 1926 1928 1930 2227 2229 2231 2286 2306 2308 2310 2314 2316 2318 2380 2384 2419 2425 2426 2427 2431 2432 2434 2436 2439 2445 2489 2492 2493 2495 2497 2498 2500 2502 2555 2606 2609 2610 2612 2614 2615 2616 2618 2620 2622 2703 2705 2706 2710")
FILENAMES=("proa" "prob" "proc" "pocket")


[ ! -z "PREFIX" ] && echo "mkvmd> Please specify prefix for output files"

if [ ! -n "$REF" ]; then
    REF="$PDB"
fi

files=("$PDB" "$TRJ" "$REF")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

cat > tcl <<EOF
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set num_frames [molinfo top get numframes]
set sel_all [atomselect top all]
mol new $REF waitfor all
set sel_ref0 [atomselect 1 "${SELREF}" frame 0]
set sel_ref [atomselect 0 "${SELREF}"]
EOF

counter=0
for region in "${REGIONS[@]}"; do
    filename="${PREFIX}-${FILENAMES[counter]}"
    rm $filename
    cat >> tcl <<EOF
set sel_rmsd0 [atomselect 1 "${region}" frame 0]
set sel_rmsd [atomselect 0 "${region}"]
set outfile [open $filename "w"]
puts \$outfile "# reference: $SELREF; calculated region: $region"
for {set i 0} {\$i<\$num_frames} {incr i} {
    \$sel_all frame \$i
    \$sel_ref frame \$i
    \$sel_rmsd frame \$i
    \$sel_all move [measure fit \$sel_ref \$sel_ref0]
    set rmsd [measure rmsd \$sel_rmsd \$sel_rmsd0 weight mass]
    puts \$outfile "\$i \t \$rmsd"
}
close \$outfile
EOF
((counter++))
done

cat >> tcl <<EOF
quit
EOF
vmd -dispdev text -e tcl

rm tcl
