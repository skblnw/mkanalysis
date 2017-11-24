set molnum [mol new ../first/initial.psf waitfor all]
mol addfile ../first/initial.pdb
mol addfile ../combine_dcd/cls2-gs3-water-protein-47ns-s50.dcd waitfor all

set nframe [molinfo $molnum get numframes]

for {set i 0} {$i < $nframe} {incr i} {
    foreach {name a b c d} {i1a 1 2 3 4 i1b 5 6 7 8 i2a 5 6 1 2 i2b 3 4 7 8} {
        set sel1 [atomselect top "segname P${a} P${b} and name CA"]
        set sel2 [atomselect top "segname P${c} P${d} and name CA"]
        $sel1 frame $i
        $sel2 frame $i
        puts "Frame $i"
        
        set res [measure contacts 10 $sel1 $sel2]
        
        set col1 [lindex $res 0]
        set col2 [lindex $res 1]
        
        set outf [open ./output/contact-$name-$i.dat "w"]
        puts "Writing results to ./output/contact-$name-$i.dat"
        for {set j 0} {$j < [llength $col1]} {incr j} {
            set idx1 [lindex $col1 $j]
            set idx2 [lindex $col2 $j]
            
            set sel3 [atomselect top "index $idx1"]
            set out1 [$sel3 get resid]
            set out3 [$sel3 get segname]
            set sel4 [atomselect top "index $idx2"]
            set out2 [$sel4 get resid]
            set out4 [$sel4 get segname]
            
            puts $outf "$out1\t$out2\t$out3\t$out4"
        }
        close $outf
    }
}
quit
