#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (contact@skblnw.com) Dec 2015
## Usage: understand, modify and bash it!
## Units: A
#########################################

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]"; exit 1; }

files=("$PDB" "$TRJ")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nFile not found!"
        exit 1
    fi
done

SELTEXT_OF_PROTEIN="backbone and resid 443 to 450"

echo "" > ${OUTPUT}.csv
cat > tcl << EOF
proc measureRgyr { nn } {

    set sel [atomselect top "${SELTEXT_OF_PROTEIN}" frame \$nn]

    return [measure rgyr \$sel]
}

# /------------------/
# /     Main Body    /
# /------------------/

# Load your structure and frames
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    
    # Define output filename
    set outf [open ${OUTPUT}.csv "a"]

    set output [measureRgyr \$nn]

    foreach element \$output {puts -nonewline \$outf [format "%.2f" \$element]}
    puts \$outf ""

    close \$outf
}
quit
EOF

vmd -dispdev text -e tcl
rm tcl
