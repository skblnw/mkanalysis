proc angle { a b } {
  # get Pi 
  global M_PI

   # Angle between two vectors 
  set cosine [expr [vecdot $a $b] / ( [veclength $a] * [veclength $b])]
  return [expr acos($cosine)*(180.0/$M_PI)]
}

set output_name [lindex $argv 0]

mol new ../cls2-combine-dcd/initial.psf
mol addfile ../cls2-combine-dcd/initial.pdb waitfor all
mol addfile ../cls2-combine-dcd/mdff-full.dcd step 2 waitfor all

set nfram [molinfo 0 get numframes]

# Script to find rotation angles between principal axes 
package require Orient 
namespace import Orient::orient 

set outf1 [open ${output_name}.dat w]

for {set j 0} {$j < $nfram } { incr j } { 
    puts "Frame $j"

    set time [expr $j*1e-15*1000*100*2*1e9]

    set output_one_line1 [list [format "%0.2f" $time]]
    foreach {k l} {1 2 3 4 5 6 7 8} {
        set sel1 [atomselect top "protein and segname P${k} and resid 1 to 210 or segname P${l} and resid 211 to 251" frame $j]
        set Isel1 [Orient::calc_principalaxes $sel1]
        set sel2 [atomselect top "protein and segname P${l} and resid 1 to 210 or segname P${k} and resid 211 to 251" frame $j]
        set Isel2 [Orient::calc_principalaxes $sel2]

        set aout [angle [lindex $Isel1 2] [lindex $Isel2 2]]
        set aout [expr 180 - $aout]
        lappend output_one_line1 [format "%0.8f" $aout]
    }

    puts $outf1 "$output_one_line1"
}
close $outf1
quit
