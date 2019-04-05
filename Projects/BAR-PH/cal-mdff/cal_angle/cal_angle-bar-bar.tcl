proc angle { a b } {
    global M_PI
    
    set cosine [expr [vecdot $a $b] / ( [veclength $a] * [veclength $b])]
    return [expr acos($cosine)*(180.0/$M_PI)]
}

set output_name [lindex $argv 0]

set molnum [mol new ../combine_dcd/cls2-4dimer.psf type psf waitfor all]
mol addfile ../combine_dcd/initial.pdb waitfor all
mol addfile ../combine_dcd/cls2-gs3-water-protein-step100-s50.dcd type dcd waitfor all

set num_frames [molinfo 0 get numframes]
set a0 [atomselect top "protein" frame 0]
set b [atomselect top "protein"]

#set outrmsd [open "result_rmsd.dat" "w"]
set outang [open output/$output_name.dat "w"]

for {set n 0} {$n < $num_frames} {incr n} {
    set time [format "%0.2f" [expr $n*1e-15*1000*100*1e9]]

    $b frame $n
    $b update
    puts "Frame $n"

#    $b move [measure fit $b $a0]
#    set rmsd [measure rmsd $b $a0]
#    puts $outrmsd "$n\t$rmsd"
    
    set ang 0
    foreach {i j k l} {1 2 3 4 5 6 7 8} {
        set selPH1 [atomselect top "segname P$i and resid 251 to 364" frame $n]
        set selPH2 [atomselect top "segname P$j and resid 251 to 364" frame $n]
        set C1 [measure center $selPH1]
        set C2 [measure center $selPH2]
        set vec1 [vecsub $C2 $C1]
        
        set selPH1 [atomselect top "segname P$k and resid 251 to 364" frame $n]
        set selPH2 [atomselect top "segname P$l and resid 251 to 364" frame $n]
        set C1 [measure center $selPH1]
        set C2 [measure center $selPH2]
        set vec2 [vecsub $C2 $C1]
        
        set ang1 [angle $vec1 $vec2]
        set ang [expr $ang + $ang1]
    }
    puts $outang "$time\t[expr $ang/2]"
}
#close $outrmsd
#close $outccc
close $outang

quit
