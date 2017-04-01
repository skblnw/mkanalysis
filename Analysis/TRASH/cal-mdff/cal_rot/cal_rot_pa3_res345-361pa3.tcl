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
mol addfile ../cls2-combine-dcd/mdff-full.dcd waitfor all

set nfram [molinfo 0 get numframes]

# Script to find rotation angles between principal axes 
package require Orient 
namespace import Orient::orient 

set selall [atomselect top "protein and resid 1 to 251"]
set outf [open $output_name.dat w]

for {set j 0} {$j < $nfram } { incr j } { 
    $selall frame $j
    $selall update
    puts "Frame $j"

    set time [expr $j*1e-15*1000*100*1e9]

#    set Ax [Orient::orient $selall [lindex $Iall 0] [lindex $Isel1 0]]
#    set Ay [Orient::orient $selall [lindex $Iall 1] [lindex $Isel1 1]]
#    set Az [Orient::orient $selall [lindex $Iall 2] [lindex $Isel1 2]]

    set output_one_line [list [format "%0.2f" $time]]
    foreach k {1 2 3 4 5 6 7 8} {
        if {$k==1 | $k==3 | $k==5 | $k==7} {
            set l [expr $k + 1]
            set selseg [atomselect top "protein and segname P${k} P${l}"]
            $selseg frame $j
            $selseg update
            set Ip [Orient::calc_principalaxes $selseg]
        }

        set sel1 [atomselect top "protein and segname P${k} and resid 345 to 361"]
        $sel1 frame $j
        $sel1 update
        set Isel1 [Orient::calc_principalaxes $sel1]
        set aout($k) [angle [lindex $Ip 2] [lindex $Isel1 2]]

        lappend output_one_line [format "%0.8f" $aout($k)]
    }

    puts $outf "$output_one_line"
}
close $outf
quit
