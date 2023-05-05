#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################


: << 'END'
# VMD/tcl Measure Script

This script is a general framework for performing various measurements using VMD/tcl. You can modify the script to include different calculation functions, selections, and input/output files.

## Usage

```
./script.sh [PDB] [TRJ] [OUTPUT]
```

## Parameters

- `[PDB]`: The input PDB file.
- `[TRJ]`: The input trajectory file.
- `[OUTPUT]`: The output file where the measurement results will be saved.

## Selections

You can modify the following selections in the script:

- `SEL1`: Selection 1.
- `SEL2`: Selection 2.
- `SEL3`: Selection 3.

## Customizing the Measurement Function

Modify the `function_template` function in the script to include the desired calculation. The function should return a value (e.g., distance, angle, etc.) that will be added to the output file.

For example, to calculate the distance between two selections, you can use the following code inside the `function_template` function:

```
set distance [measure distance $sel1 $sel2]
return $distance
```

Remember to replace the `function_template` function in the main loop with the name of your custom function.
END

# Selection for 
SEL1=""
# Selection for 
SEL2=""
# Selection for
SEL3=""

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1\n       Selection 2: $SEL2"; exit 1; }

files=("$PDB" "$TRJ")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

rm ${OUTPUT}

# Function template
function_template() {
    # Add your function code here
    # Example: calculate distance between two selections
    # set distance [measure distance \$sel1 \$sel2]
    # return \$distance
}

cat > tcl << EOF
# /------------------/
# /     Main Body    /
# /------------------/

# Load packages for calculating principle axis (if needed)
package require Orient 
namespace import Orient::orient 

mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

proc function_template {sel1 sel2} {
$(function_template)
}

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/

    set outf [open ${OUTPUT} "a"]
    # Write TIME at the very first of a line
    set out_line [format "%d" \$nn]
    # Definition of segments
    # This determines number of columns in one output file
    foreach {seg1 seg2} {1 1} {
        # Call calc funtion you like
        lappend out_line [function_template \$sel1 \$sel2]
    }
    # Write to file
    puts \$outf "\$out_line"
    # Remember to close the file
    close \$outf
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
