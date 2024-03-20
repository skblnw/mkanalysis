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

SELTEXT_OF_WATER="name OW OH2"
SELTEXT_OF_PROTEIN="noh protein and resid 340 to 375"

echo "dist6,dist7,dist8,dist9," > ${OUTPUT}.csv
cat > tcl << EOF
proc countInterfaceWater { nn } {

    set list []
    set dlist {6 7 8 9}

    foreach dd \$dlist {
        set selwater [atomselect top "${SELTEXT_OF_WATER} and same residue as water and within \$dd of resname POPC and name \"C2.*\" \"C3.*\" and within \$dd of $SELTEXT_OF_PROTEIN" frame \$nn]
        lappend list [\$selwater num]
    }

    return \$list
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

    set output [countInterfaceWater \$nn]

    foreach element \$output {puts -nonewline \$outf [format "%d," \$element]}
    puts \$outf ""

    close \$outf
}
quit
EOF

vmd -dispdev text -e tcl
rm tcl
