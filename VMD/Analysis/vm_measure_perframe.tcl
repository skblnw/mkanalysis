#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (cchan2242-c@my.cityu.edu.hk) Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

# Protocol:
# for numframes
#     for selections
#         open file for writing
#         [for segments]
#             call measure scripts
#         [end]
#     [end]
# [end]

source 

# /------------------/
# /     Main Body    /
# /------------------/

# !!!Important!!!
# Deleting existing files as we APPEND instead of trashing and opening new files
# Make sure you will delete all the existing files
# !!!Important!!!
#eval file delete [glob output/*.dat]
set OUTPUT_DIR output
exec mkdir -p $OUTPUT_DIR

# Load packages for calculating principle axis (if needed)
#package require Orient 
#namespace import Orient::orient 

# Load your structure and frames
set molnum [mol new ../combine_dcd/initial-noW.psf waitfor all]
mol addfile ../combine_dcd/initial-noW.pdb waitfor all
mol addfile ../combine_dcd/ waitfor all

set total_frame [molinfo $molnum get numframes]
for {set nn 0} {$nn < $total_frame} {incr nn} {
    set nframe [expr $nn + 0]
    puts "Frame $nframe"
    set outf [open $OUTPUT_DIR/frame$nframe.dat "a"]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Uncomment if you need PA
#    set selref [atomselect top "protein and name CA" frame $nn]
#    set Iref [Orient::calc_principalaxes $selref]
    
    foreach {sel_input1} {280} {
      set out_line [format "%d" $sel_input1]
      # This determines <number of columns in one output file>
      foreach {sel_input2} {"-19"} {
          # Call calc funtion you like
          lappend out_line [calc_dist $nn $sel_input1 $sel_input2]
      }
      
      # Write to file
      puts $outf "$out_line"
    }
    # Remember to close the file
    close $outf
}

quit
