#!/bin/bash
# KevC @ 2019 @ OSU
# Bash script to calculate distance between two residues

PDB="$1"
TRJ="$2"
OPT="$3"
SEL1="$4"
SEL2="$5"
[ $# -ne 5 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [-i|-r] [Selection 1] [Selection 2]"; echo "mkvmd> By default, the selection is 'protein and resid ${sel_input1} and name CA'"; exit 1; }
case $OPT in
    -r)
    SELTEXT1="protein and resid $SEL1 and name CA"
    SELTEXT2="protein and resid $SEL2 and name CA"
    ;;
    -i)
    SELTEXT1="index $SEL1"
    SELTEXT2="index $SEL2"
    ;;
    *)    # unknown option
    echo -e "Unknown option $OPT, must be either -i or -r"
    exit 0
    ;;
esac


if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

echo "" > vm_cal-dist.tcl
cat > vm_cal-dist.tcl << EOF

proc calc_dist {nn} {
set sel1 [atomselect top "$SELTEXT1" frame \$nn]
set sel2 [atomselect top "$SELTEXT2" frame \$nn]
set dist [vecdist [measure center \$sel1 weight mass] [measure center \$sel2 weight mass]]
set res [format "%.2f" \$dist]
\$sel1 delete
\$sel2 delete
return \$res
}
set molnum [mol new $PDB waitfor all]
mol addfile $TRJ waitfor all
set total_frame [molinfo \$molnum get numframes]
set outf [open out.dat w]
for {set nn 0} {\$nn < \$total_frame} {incr nn} {puts \$outf "[expr \$nn + 1] [calc_dist \$nn]"}
close \$outf
quit

EOF

vmd -dispdev text -e vm_cal-dist.tcl
