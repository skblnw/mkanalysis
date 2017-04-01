set filename [lindex $argv 0]
set outname [lindex $argv 1]


set molnum [mol new ../combine_dcd/cls1-4dimer.psf type psf waitfor all]
mol addfile ../combine_dcd/cls1-gs3-water-protein-8ns-s50.dcd type dcd waitfor all

puts "M) Number of frames read: "
set num_frames [molinfo $molnum get numframes]
set outf [open output/$filename-$outname.dat "w"]

for {set i 0} {$i < $num_frames} {incr i} {
    set total1 0
    set total2 0
    foreach {a b c d} {1 2 5 6 7 8 3 4} {
        set seg1 "P$a P$b"
        set seg2 "P$c P$d"
    
        set sel1 [atomselect top "segname $seg1 and name CA"]
        set sel2 [atomselect top "segname $seg2 and name CA"]
        $sel1 frame $i
        $sel2 frame $i
        
        set res [measure contacts 10 $sel1 $sel2]
        
        set col1 [lindex $res 0]
        set col2 [lindex $res 1]
        set num1 [llength $col1]
        set num2 [llength $col2]

        set total1 [expr $total1 + $num1]
        set total2 [expr $total2 + $num2]

        if [expr $i == $num_frames-1] {
            set out_full [open output/$filename-$outname.full "a"]
            for {set j 0} {$j < [llength $col1]} {incr j} {
          
              set idx1 [lindex $col1 $j]
              set idx2 [lindex $col2 $j]
          
              set sel3 [atomselect top "index $idx1"]
              set out1 [$sel3 get resid]
              set out3 [$sel3 get segname]
              set sel4 [atomselect top "index $idx2"]
              set out2 [$sel4 get resid]
              set out4 [$sel4 get segname]
          
              puts $out_full "$out1\t$out2\t$out3\t$out4"
            }
            close $out_full
        }

    }
    
    if {$total1 == $total2} {
        puts $outf "$i\t$total1"
    } else {
        puts $outf "ERROR: [Step $i] Columns don't match"
    }

}

close $outf
quit
