#!/bin/bash
# KevC @ 2019 @ OSU
# Bash script to calculate distance between two residues

[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [PREFIX] [-i|-r] [Selection 1] [Selection 2]"; exit 1; }

PDB="$1"
TRJ="$2"
PREFIX="$3"
OPT="$4"
SEL1="$5"
SEL2="$6"

files=("$PDB" "$TRJ")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

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
set molnum [mol new $PDB waitfor all]
mol addfile $TRJ waitfor all
set total_frame [molinfo \$molnum get numframes]
set outf [open "${PREFIX}" w]
EOF

case $OPT in
    -r)
        SELTEXT1="protein and resid $SEL1 and name CA"
        SELTEXT2="protein and resid $SEL2 and name CA"
        cat >> tcl <<EOF
for {set nn 0} {\$nn < \$total_frame} {incr nn} {puts \$outf "\$nn [calc_dist \$nn \"$SELTEXT1\" \"$SELTEXT2]\""}
close \$outf
quit
EOF
    ;;
    -i)
        cat >> tcl <<EOF
puts \$outf "# option: $OPT; selection 1: $SEL1; selection 2: $SEL2"
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
