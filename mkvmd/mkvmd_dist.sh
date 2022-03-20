#!/bin/bash
# KevC @ 2019 @ OSU
# Bash script to calculate distance between two residues
source ~/.zshrc

PDB="$1"
TRJ="$2"
OUTPUT="$3"
OPT="$4"
SEL1="$5"
SEL2="$6"
[ $# -ne 6 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [-i|-r] [Selection 1] [Selection 2]"; echo "mkvmd> By default, the selection is 'protein and resid ${sel_input1} and name CA'"; exit 1; }


if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

cat > tcl << EOF
proc calc_dist {nn sel_input1 sel_input2} {
    set sel1 [atomselect top "\$sel_input1" frame \$nn]
    set sel2 [atomselect top "\$sel_input2" frame \$nn]
    set dist [vecdist [measure center \$sel1 weight mass] [measure center \$sel2 weight mass]]
    set res [format "%.2f" \$dist]
    \$sel1 delete
    \$sel2 delete
    return \$res
}
EOF

case $OPT in
    -r)
        SELTEXT1="protein and resid $SEL1 and name CA"
        SELTEXT2="protein and resid $SEL2 and name CA"
        cat >> tcl <<EOF
set molnum [mol new $PDB waitfor all]
mol addfile $TRJ waitfor all
set total_frame [molinfo \$molnum get numframes]
set outf [open dist w]
for {set nn 0} {\$nn < \$total_frame} {incr nn} {puts \$outf "[expr \$nn + 1] [calc_dist \$nn]"}
close \$outf
quit
EOF
    ;;
    -i)
        cat >> tcl <<EOF
set molnum [mol new $PDB waitfor all]
mol addfile $TRJ waitfor all
set total_frame [molinfo \$molnum get numframes]
set outf [open "$OUTPUT" w]
for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set out_line [format "%d" \$nn]
    foreach {idx1} {$SEL1} {idx2} {$SEL2} {
        lappend out_line [calc_dist \$nn "index \$idx1" "index \$idx2"]
    }
    puts \$outf "\$out_line"
}
close \$outf
quit
EOF
    ;;
    *)    # unknown option
        echo -e "Unknown option $OPT, must be either -i or -r"
        exit 0
    ;;
esac
vmd -dispdev text -e tcl
rm tcl
