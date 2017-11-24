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

source proc_findContacts.tcl

# /------------------/
# /     Main Body    /
# /------------------/

# !!!Important!!!
# Deleting existing files as we APPEND instead of trashing and opening new files
# Make sure you will delete all the existing files
# !!!Important!!!
#eval file delete [glob output/*.dat]
set OUTPUT_DIR output_com_95100_12_frame
exec mkdir -p $OUTPUT_DIR

# Load packages for calculating principle axis (if needed)
#package require Orient 
#namespace import Orient::orient 

# Load your structure and frames
set molnum [mol new ../combine_dcd/initial-noW.psf waitfor all]
mol addfile ../combine_dcd/initial-noW.pdb waitfor all
mol addfile ../combine_dcd/npt-pf50ps-13.dcd waitfor all

set total_frame [molinfo $molnum get numframes]
for {set nn 0} {$nn < $total_frame} {incr nn} {
    set nframe [expr $nn + 0]
    puts "Frame $nframe"

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Uncomment if you need PA
#    set selref [atomselect top "protein and name CA" frame $nn]
#    set Iref [Orient::calc_principalaxes $selref]

    # Reset index of selections
    set nsel 0
    
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
    # This determines <number of output files>
    foreach {sel_input1} {1} {sel_input2} {1} {
      # Output file name
      # Increase index of selections by 1
      incr nsel
      set outf [open $OUTPUT_DIR/sel$nsel.dat "a"]
      # Write TIME at the very first of a line
      set out_line [format "%d" $nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
    # This determines <number of columns in one output file>
      foreach {seg1 seg2} {1 1} {
          # Call calc funtion you like
          lappend out_line [findContacts $nn 12]
      }
      # Write to file
      puts $outf "$out_line"
      # Remember to close the file
      close $outf
    }
}

quit
