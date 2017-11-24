proc angle { a b } {
    global M_PI
    
    set cosine [expr [vecdot $a $b] / ( [veclength $a] * [veclength $b])]
    return [expr acos($cosine)*(180.0/$M_PI)]
}

set output_name [lindex $argv 0]

set molnum [mol new ../combine_dcd/initial.psf type psf waitfor all]
mol addfile ../combine_dcd/initial.pdb waitfor all
mol addfile ../combine_dcd/mdff-s50.dcd type dcd waitfor all

set outf [open $output_name.dat "w"]

set num_frames [molinfo 0 get numframes]
for {set n 0} {$n < $num_frames} {incr n} {
    set time [expr $n*1e-15*1000*100*1e9]

    set output_one_line1 [list [format "%0.2f" $time]]
    foreach {i j k} {19 20 7 7 8 3 1 2 5 5 6 17 23 24 15 15 16 11 9 10 13 13 14 21} {
        set selPH1 [atomselect top "segname P${i} and resid 251 to 364" frame $n]
        set selPH2 [atomselect top "segname P${j} and resid 251 to 364" frame $n]
        set selPH3 [atomselect top "segname P${k} and resid 251 to 364" frame $n]
        set C1 [measure center $selPH1]
        set C2 [measure center $selPH2]
        set C3 [measure center $selPH3]
        set line1 [vecsub $C2 $C1]
        set line2 [vecsub $C3 $C1]
        
        set aout [angle $line1 $line2]
        lappend output_one_line1 [format "%0.4f" $aout]
    }
    puts $outf "$output_one_line1"
}
close $outf

quit