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

source proc_sel2norm.tcl
source proc_calc-2vec.tcl

# /------------------/
# /     Main Body    /
# /------------------/

# !!!Important!!!
# Deleting existing files as we APPEND instead of trashing and opening new files
# Make sure you will delete all the existing files
# !!!Important!!!
#eval file delete [glob output/*.dat]
set OUTPUT_DIR output_pa
exec mkdir -p $OUTPUT_DIR

# Load packages for calculating principle axis (if needed)
package require Orient 
namespace import Orient::orient 

# Load your structure and frames
set molnum [mol new ../combine_dcd/initial-noW.psf waitfor all]
mol addfile ../combine_dcd/initial-noW.pdb waitfor all
mol addfile ../combine_dcd/npt-pf1000ps-13.dcd waitfor all

set total_frame [molinfo $molnum get numframes]
for {set nn 0} {$nn < $total_frame} {incr nn} {
    set nframe [expr $nn + 0]
    puts "Frame $nframe"

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Uncomment if you need PA
    set selref [atomselect top "protein and resid 1 to 250 and name CA" frame $nn]
    set Iref [Orient::calc_principalaxes $selref]
    set Iref1 [lindex $Iref 0]
    set Iref2 [lindex $Iref 1]
    set Iref3 [lindex $Iref 2]
    if { [angle $Iref1 {0 1 0}] > 90 } {
      set Iref1 [vecscale $Iref1 -1]
    }
    if { [angle $Iref2 {0 0 1}] > 90 } {
      set Iref2 [vecscale $Iref2 -1]
    }
    if { [angle $Iref3 {-1 0 0}] > 90 } {
      set Iref3 [vecscale $Iref3 -1]
    }

    # Reset index of selections
    set nsel 0
    
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
    # This determines <number of output files>
    foreach {sel_input1} {"277 to 283" "277 to 283"} {sel_input2} {"321 to 324" "300 to 302"} {sel_input3} {"325 to 327" "303 to 305"} {
      # Output file name
      # Increase index of selections by 1
      incr nsel
      set outf1 [open $OUTPUT_DIR/sel$nsel-BAR1.dat "a"]
      set outf2 [open $OUTPUT_DIR/sel$nsel-BAR2.dat "a"]
      set outf3 [open $OUTPUT_DIR/sel$nsel-BAR3.dat "a"]
      # Write TIME at the very first of a line
      set out_line1 [format "%d" $nframe]
      set out_line2 [format "%d" $nframe]
      set out_line3 [format "%d" $nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
    # This determines <number of columns in one output file>
      foreach {seg} {1 2} {
          # Call calc funtion you like
          set vec [sel2norm $nn $seg $seg $seg $sel_input1 $sel_input2 $sel_input3]

          lappend out_line1 [calc_2vec $nn $vec $Iref1]
          lappend out_line2 [calc_2vec $nn $vec $Iref2]
          lappend out_line3 [calc_2vec $nn $vec $Iref3]
      }
      # Write to file
      puts $outf1 "$out_line1"
      puts $outf2 "$out_line2"
      puts $outf3 "$out_line3"
      # Remember to close the file
      close $outf1
      close $outf2
      close $outf3
    }
}

quit
