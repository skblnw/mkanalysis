#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

SEL_PEPTIDE="segname PROC"
SEL_HLA="segname PROA and resid 1 to 180"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1\n       Selection 2: $SEL2\n       Selection combined: $SEL_COMBINE"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

 # and not name CA C N O HN HA

echo "pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9," > ${OUTPUT}_mainnp.csv
echo "pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9," > ${OUTPUT}_mainp.csv
echo "pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9," > ${OUTPUT}_sidenp.csv
echo "pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9," > ${OUTPUT}_sidep.csv
cat > tcl << EOF
proc countIntermolContact { nn sel_input } {
    set sel [atomselect top "$SEL_PEPTIDE"]
    set reslist [lsort -unique -integer [\$sel get residue]]
    if {\$nn == 0} {
        puts "Selected [llength \$reslist] residues"
    }

    set output []
    set dlist {3.5}
    foreach ii \$reslist {
        set tmplist 0
        foreach dd \$dlist {
            set selheavy [atomselect top "\$sel_input and not residue \$ii and within \$dd of {residue \$ii}" frame \$nn]
            incr tmplist [\$selheavy num]
        }
        lappend output [expr \$tmplist / [llength \$dlist]]
    }
    
    set total 0
    foreach nxt \$output { incr total \$nxt }
    puts "Total number of nearby atoms: \$total"
    return \$output
}

# /------------------/
# /     Main Body    /
# /------------------/

mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/

    # Output file name
    set outf1 [open ${OUTPUT}_mainnp.csv "a"]
    set outf2 [open ${OUTPUT}_mainp.csv "a"]
    set outf3 [open ${OUTPUT}_sidenp.csv "a"]
    set outf4 [open ${OUTPUT}_sidep.csv "a"]

    # Call calc funtion you like
    set output1 [countIntermolContact \$nn "$SEL_HLA and backbone and name \"C.*\""]
    set output2 [countIntermolContact \$nn "$SEL_HLA and backbone and name \"O.*\" \"N.*\""]
    set output3 [countIntermolContact \$nn "$SEL_HLA and sidechain and name \"C.*\" \"S.*\""]
    set output4 [countIntermolContact \$nn "$SEL_HLA and sidechain and name \"O.*\" \"N.*\""]

    # Write to file
    foreach element \$output1 {puts -nonewline \$outf1 "\$element,"}
    puts \$outf1 ""
    foreach element \$output2 {puts -nonewline \$outf2 "\$element,"}
    puts \$outf2 ""
    foreach element \$output3 {puts -nonewline \$outf3 "\$element,"}
    puts \$outf3 ""
    foreach element \$output4 {puts -nonewline \$outf4 "\$element,"}
    puts \$outf4 ""

    # Remember to close the file
    close \$outf1
    close \$outf2
    close \$outf3
    close \$outf4
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
