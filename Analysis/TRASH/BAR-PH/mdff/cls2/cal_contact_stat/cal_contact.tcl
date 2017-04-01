set filename [lindex $argv 0]
set seltext1 [lrange $argv 1 3]
set seltext2 [lrange $argv 4 6]
set outname [lindex $argv 7]


set molnum [mol new ../combine_dcd/cls2-4dimer.psf type psf waitfor all]
mol addfile ../combine_dcd/cls2-gs3-water-protein-11ns-s50.dcd type dcd waitfor all

puts "M) Number of frames read: "
set num_frames [molinfo $molnum get numframes]
set outf [open $filename-$outname.dat "w"]

for {set i 0} {$i < $num_frames} {incr i} {
    set total1 0
    set total2 0
    foreach {a b c d} {5 6 1 2 3 4 7 8} {
        set seg1 "P$a P$b"
        set seg2 "P$c P$d"
    
        set sel1 [atomselect top "segname $seg1 and resid $seltext1 and name CA"]
        set sel2 [atomselect top "segname $seg2 and resid $seltext2 and name CA"]
        $sel1 frame $i
        $sel2 frame $i
        
        set res [measure contacts 10 $sel1 $sel2]
        
        set col1 [lindex $res 0]
        set col2 [lindex $res 1]
        set num1 [llength $col1]
        set num2 [llength $col2]

        set total1 [expr $total1 + $num1]
        set total2 [expr $total2 + $num2]

    }
    
    if {$total1 == $total2} {
        puts $outf "$i\t$total1"
    } else {
        puts $outf "ERROR: [Step $i] Columns don't match"
    }

}

close $outf
quit
