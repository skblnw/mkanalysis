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

SELTEXT1="segname PROC"
SELTEXT2="name OW"

rm $OUTPUT
cat > tcl << EOF
proc lcompare { list1 list2 } {
    set outlist {}
    foreach ii \$list1 {
        if {[lsearch -exact \$list2 \$ii] != -1} {
            lappend outlist \$ii
        }
    }
    if {[llength \$outlist] != 0} {
        return \$outlist
    } else {
        return -1
    }
}

proc countHetwork {nn donor acceptor} {
    set outlist {}

    set selsc [atomselect top "${SELTEXT1} and resid \$donor \$acceptor and not backbone" frame \$nn]
    set selw1 [atomselect top "${SELTEXT2} and within 3.5 of ${SELTEXT1} and resid \$donor and not backbone" frame \$nn]
    set selw2 [atomselect top "${SELTEXT2} and within 3.5 of ${SELTEXT1} and resid \$acceptor and not backbone" frame \$nn]

    set list1 [measure hbonds 3.5 30 \$selsc \$selw1]
    set list2 [measure hbonds 3.5 30 \$selw1 \$selw2]
    set list3 [measure hbonds 3.5 30 \$selw2 \$selsc]

    set list4 [measure hbonds 3.5 30 \$selsc]

    if {[llength [lindex \$list4 0]] != 0} {
        lappend outlist 1
    } else {
        lappend outlist 0
    }

    if {[lcompare [lindex \$list1 1] [lindex \$list3 0]] != -1} {
        lappend outlist 1
    } else {
        lappend outlist 0
    }

    foreach ii [lindex \$list2 0] jj [lindex \$list2 1] kk [lindex \$list2 2] {
        set counter 0
        if {[lsearch -exact [lcompare [lindex \$list1 1] [lindex \$list2 0]] \$ii] != -1} {incr counter}
        if {[lsearch -exact [lcompare [lindex \$list2 1] [lindex \$list3 0]] \$jj] != -1} {incr counter}
        if { \$counter != 2} {
            continue
        } else {
            lappend outlist 1
            break
        }
    }

    if {[llength \$outlist] == 2} {
        lappend outlist 0
        return \$outlist
    } else {
        return \$outlist
    }
}

# /------------------/
# /     Main Body    /
# /------------------/

# /------------------------------------------------/
# / Deleting existing files as we APPEND instead of trashing and opening new files
# /------------------------------------------------/
# eval file delete [glob output/*.dat]
# set OUTPUT_DIR [exec date +%Y%m%d%H%M%S]
# exec mkdir -p \$OUTPUT_DIR

# Load packages for calculating principle axis if needed
# package require Orient 
# namespace import Orient::orient 

# Load your structure and frames
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]
    set out_line [format "%d" \$nframe]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Uncomment these two lines if you need PA
    # set selref [atomselect top "protein and name CA" frame \$nn]
    # set Iref [Orient::calc_principalaxes \$selref]
    # /------------------------------------------------/
    # /                 Atom Selections
    # / You may use "...", e.g. "1 to 10", instead of one integer
    # / This determines <number of output files>
    # /------------------------------------------------/
    foreach {sel_input1} {493} {sel_input2} {484} {
        set outf [open $OUTPUT "a"]
        # Write to file
        foreach ii [countHetwork \$nn \$sel_input1 \$sel_input2] {lappend out_line \$ii}
        puts \$outf "\$out_line"
        # Remember to close the file
        close \$outf
    }
}

quit
EOF

vmd -dispdev text -e tcl
rm -f tcl
