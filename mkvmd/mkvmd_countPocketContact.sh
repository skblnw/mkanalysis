#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

# Selection for peptide
SEL_PEPTIDE="segname PROC"
SEL_HLA="segname PROA"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1\n       Selection 2: $SEL2"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

rm ${OUTPUT}_*
cat > tcl << EOF
proc countPocketContact { frame_number count_dict protein_residues_selection peptide_residues } {
    set peptide_selection [atomselect top "$SEL_PEPTIDE"]

    foreach peptide_residue \$peptide_residues {
        set heavy_atom_selection [atomselect top "noh \$protein_residues_selection and within 4 of $SEL_PEPTIDE and resid \$peptide_residue" frame \$frame_number]
        if { [\$heavy_atom_selection num] != 0 } {
            dict incr count_dict \$peptide_residue
        }
    }
    return [dict values \$count_dict]
}

# /------------------/
# /     Main Body    /
# /------------------/

# Load packages for calculating principle axis (if needed)
package require Orient 
namespace import Orient::orient 

# Load trajectory
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

# Define peptide_residues (generally 1-9)
set peptide_selection [atomselect top "$SEL_PEPTIDE"]
set peptide_residues [lsort -unique -integer [\$peptide_selection get resid]]
puts "Selected [llength [lsort -unique -integer [\$peptide_selection get residue]]] residues"

# Define protein_residues (list from zql)
set protein_residues_selections [list]
foreach ii { 5 7 9 24 33 34 45 55 59 60 62 63 65 66 67 69 70 72 73 74 76 77 80 81 83 84 95 97 99 114 116 118 123 126 133 143 146 147 150 152 155 156 158 159 160 163 167 170 171 } {
    lappend protein_residues_selections "$SEL_HLA and resid \$ii"
}

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {

    foreach protein_residues_selection \$protein_residues_selections {

        set total 0

        # Define output filename based on the protein residue selection
        set output_filename [string map {" " "_" "=" ""} \$protein_residues_selection]
        set output_file "${OUTPUT}_\${output_filename}"

        # Open file for writing
        set outf [open \$output_file a]

        # Initialize count_dict for each protein residue selection
        set count_dict [dict create]
        foreach residue \$peptide_residues {
            dict set count_dict \$residue 0
        }

        set output [countPocketContact \$nn \$count_dict \$protein_residues_selection \$peptide_residues]
        puts \$outf "\$output"

        foreach nxt \$output { incr total \$nxt }
        puts -nonewline "\$total "

        # Remember to close the file
        close \$outf

    }

    puts ""

}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
