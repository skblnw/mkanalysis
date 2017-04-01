proc angle { a b } {
  # get Pi 
  global M_PI

   # Angle between two vectors 
  set cosine [expr [vecdot $a $b] / ( [veclength $a] * [veclength $b])]
  return [expr acos($cosine)*(180.0/$M_PI)]
}

set output_name [lindex $argv 0]

set molnum [mol new ../combine_dcd/cls2-4dimer.psf type psf waitfor all]
mol addfile ../combine_dcd/initial.pdb waitfor all
mol addfile ../combine_dcd/cls2-gs3-water-protein-47ns-s50.dcd type dcd waitfor all

set nfram [molinfo 0 get numframes]

# Script to find rotation angles between principal axes 
package require Orient 
namespace import Orient::orient 

set outf1 [open output/${output_name}1.dat w]
set outf2 [open output/${output_name}2.dat w]
set outf3 [open output/${output_name}3.dat w]

for {set j 0} {$j < $nfram } { incr j } { 
    puts "Frame $j"

    set time [expr $j*1e-15*1000*100*1e9]

    set output_one_line1 [list [format "%0.2f" $time]]
    set output_one_line2 [list [format "%0.2f" $time]]
    set output_one_line3 [list [format "%0.2f" $time]]
    foreach k {1 2 3 4 5 6 7 8} {
        if {$k==1 | $k==3 | $k==5 | $k==7} {
            set l [expr $k + 1]
            set selseg [atomselect top "protein and segname P${k} P${l} and resid 1 to 251 and name CA"]
            $selseg frame $j
            $selseg update
            set Ip [Orient::calc_principalaxes $selseg]
        }

        set sel1 [atomselect top "protein and segname P${k} and resid 345 to 361 and name CA"]
        $sel1 frame $j
        $sel1 update
        set Isel1 [Orient::calc_principalaxes $sel1]

        if { [angle [lindex $Ip 0] {0 0 1}] > 90 } {
            set vec1 [vecinvert [lindex $Ip 0]]
        } else {
            set vec1 [lindex $Ip 0]
        }

        if { [angle [lindex $Ip 1] {0 1 0}] > 90 } {
            set vec2 [vecinvert [lindex $Ip 1]]
        } else {
            set vec2 [lindex $Ip 1]
        }

        if { [angle [lindex $Ip 2] {1 0 0}] > 90 } {
            set vec3 [vecinvert [lindex $Ip 2]]
        } else {
            set vec3 [lindex $Ip 2]
        }

        if { [angle [lindex $Isel1 2] {0 0 1}] > 90 } {
            set vec4 [vecinvert [lindex $Isel1 2]]
        } else {
            set vec4 [lindex $Isel1 2]
        }

        set aout($k) [angle $vec1 $vec4]
        lappend output_one_line1 [format "%0.8f" $aout($k)]

        set aout($k) [angle $vec2 $vec4]
        lappend output_one_line2 [format "%0.8f" $aout($k)]

        set aout($k) [angle $vec3 $vec4]
        lappend output_one_line3 [format "%0.8f" $aout($k)]
    }

    puts $outf1 "$output_one_line1"
    puts $outf2 "$output_one_line2"
    puts $outf3 "$output_one_line3"
}
close $outf1
close $outf2
close $outf3
quit
