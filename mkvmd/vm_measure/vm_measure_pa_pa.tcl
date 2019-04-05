#########################################
## Description: Modify from vm_measure.tcl. It would be a bit complicated to accurately measure rotational angle between PH and BAR domains inside a protein lattice so I specifically save this.
## Author: Kev Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

source proc_calc-2vec.tcl

# /------------------/
# /     Main Body    /
# /------------------/

# !!!Important!!!
# Deleting existing files as we APPEND instead of trashing and opening new files
# Make sure you will delete all the existing files
# !!!Important!!!

# Load packages for calculating principle axis (if needed)
package require Orient 
namespace import Orient::orient

# Load your structure and frames
set molnum [mol new ../combine_dcd/initial.psf waitfor all]
mol addfile ../combine_dcd/initial.pdb waitfor all
mol addfile ../combine_dcd/mdff-cls2-pf100ps.dcd waitfor all

set total_frame [molinfo $molnum get numframes]
for {set nn 0} {$nn < $total_frame} {incr nn} {
    puts "Frame $nn"

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Reset index of selections
    set nsel 0
    # Definition of selections
    # You may use, e.g. "1 to 10", instead of one integer
    foreach {sel_input1} {"1 to 251"} {sel_input2} {"345 to 361"} {
      # Output file name
      set outf1 [open output3/ang-rot-PH-BAR-1.dat "a"]
      set outf2 [open output3/ang-rot-PH-BAR-2.dat "a"]
      set outf3 [open output3/ang-rot-PH-BAR-3.dat "a"]
      # Write TIME at the very first of a line
      set out_line1 [format "%4.2f" $nn]
      set out_line2 [format "%4.2f" $nn]
      set out_line3 [format "%4.2f" $nn]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
      foreach seg {1 2 3 4 5 6 7 8} {
          if { [expr $seg % 2] != 0 } {
            set jj [expr $seg + 1]
            set selbar [atomselect top "protein and segname P${seg} P${jj} and resid ${sel_input1} and name CA" frame $nn]
            set Ib [Orient::calc_principalaxes $selbar]; puts ""
            $selbar delete
            # You may notice that we used BAR domains of both monomers to calculate the PA of BAR
            set Ib1 [lindex $Ib 0]
            set Ib2 [lindex $Ib 1]
            set Ib3 [lindex $Ib 2]
            if { [angle $Ib1 {0 0 -1}] > 90 } {
              set Ib1 [vecscale $Ib1 -1]
            }
            if { [angle $Ib2 {0 -1 0}] > 90 } {
              set Ib2 [vecscale $Ib2 -1]
            }
            if { [angle $Ib3 {1 0 0}] > 90 } {
              set Ib3 [vecscale $Ib3 -1]
            }
            # Split PA into three axes and align with the VMD axes.
            # [Important] Note that here we usually define the "VMD axes" by looking at segment P1.
            # Yes orientation of protein inside the lattice does cause difference in defining the rotational angle but will be solved later solely by inverting (or not) the PA of PH domains.
          }
          set selph [atomselect top "protein and segname P${seg} and resid ${sel_input2} and name CA" frame $nn]
          set Ip [Orient::calc_principalaxes $selph]; puts ""
          $selph delete
          # Calculating PA of PH domains
          set Ip1 [lindex $Ip 0]
          set Ip2 [lindex $Ip 1]
          set Ip3 [lindex $Ip 2]
          if { [angle $Ip3 {0 0 -1}] > 90 } {
            set Ip3 [vecscale $Ip3 -1]
          }

          # Again align with the VMD axis by looking at segment P1
          # [Important] Here comes the tricky part! Get a paper and open VMD, see whether you have to invert or not for each alpha angle.
          set alpha2 [list 2 3 6 7]
          lappend out_line1 [calc_2vec $Ib1 $Ip3]
          if [lsearch $alpha2 $seg]!=-1 {
              lappend out_line2 [calc_2vec $Ib2 [vecinvert $Ip3]]
          } else {
              lappend out_line2 [calc_2vec $Ib2 $Ip3]
          }
          lappend out_line3 [calc_2vec $Ib3 $Ip3]
      }
      # Write to file
      puts $outf1 "$out_line1"
      puts $outf2 "$out_line2"
      puts $outf3 "$out_line3"
      # Remember to close the file
      close $outf1
      close $outf2
      close $outf3
      incr nsel
    }
}

quit
