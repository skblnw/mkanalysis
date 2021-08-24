proc angle { a b } {
  # get Pi 
  global M_PI

   # Angle between two vectors 
  set cosine [expr [vecdot $a $b] / ( [veclength $a] * [veclength $b])]
  return [expr acos($cosine)*(180.0/$M_PI)]
}

set output_name [lindex $argv 0]

mol new ../combine_dcd/initial.psf
mol addfile ../combine_dcd/initial.pdb
mol addfile ../combine_dcd/mdff-s50.dcd waitfor all

set nfram [molinfo 0 get numframes]

# Script to find rotation angles between principal axes 
package require Orient 
namespace import Orient::orient 

set outf1 [open ${output_name}.dat w]

for {set j 0} {$j < $nfram } { incr j } { 
    puts "Frame $j"

    set time [expr $j*1e-15*1000*50*1e9]

    set output_one_line1 [list [format "%0.2f" $time]]
    foreach {k l m n} {1 2 3 4 5 6 7 8 17 18 19 20 9 10 11 12 13 14 15 16 21 22 23 24} {
        set sel1 [atomselect top "protein and segname P${k} P${l} and resid 1 to 210 and name CA" frame $j]
        set Isel1 [Orient::calc_principalaxes $sel1]
        set sel2 [atomselect top "protein and segname P${m} P${n} and resid 1 to 210 and name CA" frame $j]
        set Isel2 [Orient::calc_principalaxes $sel2]

        set aout [angle [lindex $Isel1 2] [lindex $Isel2 2]]
        if {$aout < 90} {
            set aout [expr 180 - $aout]
        }
        lappend output_one_line1 [format "%0.8f" $aout]
    }

    puts $outf1 "$output_one_line1"
}
close $outf1
quit