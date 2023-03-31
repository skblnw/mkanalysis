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
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

# SELTEXT="segname PROA and resid 56 to 89 137 to 177 and resname ARG LYS GLU ASP"
SELTEXT="segname PROC"
SELTEXT_OF_WATER="name OW OH2"

rm ${OUTPUT}_w ${OUTPUT}_hb
cat > tcl << EOF
proc countWater { nn reslist } {
    set list_water []
    set list_hbonds []
    set dlist {3 3.5}
    foreach ii \$reslist {
        set count 0
        foreach dd \$dlist {
            set selwater [atomselect top "${SELTEXT_OF_WATER} and within \$dd of residue \$ii" frame \$nn]
            incr count [\$selwater num]
        }
        lappend list_water [expr \$count / [llength \$dlist]]

        set count 0
        foreach dd \$dlist {
            set sel1 [atomselect top "residue \$ii and sidechain" frame \$nn]
            set sel2 [atomselect top "${SELTEXT_OF_WATER} and within \$dd of residue \$ii" frame \$nn]
            set count1 [llength [lindex [measure hbonds \$dd 20 \$sel1 \$sel2] 0]]
            set count2 [llength [lindex [measure hbonds \$dd 20 \$sel2 \$sel1] 0]]
            incr count [expr \$count1 + \$count2]
        }
        lappend list_hbonds [expr \$count / [llength \$dlist]]
    }
    
    set total1 0
    set total2 0
    foreach nxt \$list_water { incr total1 \$nxt }
    foreach nxt \$list_hbonds { incr total2 \$nxt }
    puts "Total number of water: \$total1. Total number of h-bonds: \$total2."
    return [list \$list_water \$list_hbonds]
}

# /------------------/
# /     Main Body    /
# /------------------/

# Load your structure and frames
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /-------------------------------------------------/
    # /     Where you really have to use your brain     /
    # /-------------------------------------------------/
    # /--------------------------------------------------------------/
    # /                       Atom Selections                        /
    # /   You may use "...", e.g. "1 to 10", instead of one integer  /
    # /          This determines <number of output columns>            /
    # /--------------------------------------------------------------/
    set outf1 [open ${OUTPUT}_w "a"]
    set outf2 [open ${OUTPUT}_hb "a"]

    set out_line1 [format "%d" \$nframe]
    set out_line2 [format "%d" \$nframe]

    set sel [atomselect top "${SELTEXT}"]
    set reslist [lsort -unique -integer [\$sel get residue]]
    if { \$nn == 0 } {
        puts "Selected [llength \$reslist] residues"
        set out_line0 "nframe "
        foreach ii \$reslist {
            set sel [atomselect top "residue \$ii and name CA"]
            lappend out_line0 "[\$sel get resname][\$sel get resid]"
        }
        puts \$outf1 "\$out_line0"
        puts \$outf2 "\$out_line0"
    }

    set output [countWater \$nn \$reslist]
    foreach ii [lindex \$output 0] {lappend out_line1 \$ii}
    foreach ii [lindex \$output 1] {lappend out_line2 \$ii}

    puts \$outf1 "\$out_line1"
    puts \$outf2 "\$out_line2"
    # Remember to close the file
    close \$outf1
    close \$outf2
}
quit
EOF

vmd -dispdev text -e tcl
rm -f tcl
