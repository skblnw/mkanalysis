#!/bin/bash
# KevC @ 2019 @ OSU
# Further Improved Bash script to calculate distance between two residues

if [ "$#" -lt 3 ]; then
    echo "mkvmd> Usage: $0 [PDB] [TRJ] [PREFIX] [INDEX_FILE]"
    exit 1
fi

PDB="$1"
TRJ="$2"
PREFIX="$3"
INDEX_FILE="$4"

# Check for existence of input files
files=("$PDB" "$TRJ" "$INDEX_FILE")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nFile not found!"
        exit 1
    fi
done

# TCL proc and common setup
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
EOF

# Process each line in the INDEX_FILE
index_count=0
while read -r IDX1 IDX2; do
    # Skip lines starting with '#' or ';'
    [[ "$IDX1" == \#* ]] || [[ "$IDX1" == \;* ]] && continue

    let index_count=index_count+1
    OUTFILE="${PREFIX}-d${index_count}"
    cat >> tcl << EOF
set outf [open "$OUTFILE" w]
puts \$outf "# Index 1: $IDX1; Index 2: $IDX2"
for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    puts \$outf "\$nn [calc_dist \$nn "index $IDX1" "index $IDX2"]"
}
close \$outf
EOF
done < "$INDEX_FILE"

cat >> tcl << EOF
quit
EOF

vmd -dispdev text -e tcl
rm tcl
