# Find salt bridges between "seltext" and "NOT seltext" within ondist
# Simply changing "seg" to two atom indexes makes the script compute atom pairs
proc findSaltBridges { nn seg sel_input1 sel_input2 ondist} {
    set sel1 [atomselect top "(protein and acidic and oxygen and not backbone) and segname D$seg and resid $sel_input1" frame $nn]
    set sel2 [atomselect top "(protein and basic and nitrogen and not backbone) and not segname D$seg and resid $sel_input2" frame $nn]
    set tmpList [measure contacts $ondist $sel1 $sel2]
    $sel1 delete
    $sel2 delete
    
    foreach i [lindex $tmpList 0] j [lindex $tmpList 1] {
        set potPairs($i,$j) 1
    }

    foreach pair [array names potPairs] {
    
        foreach { ac ba } [split $pair ,] break

        set refAc [atomselect top "same residue as index $ac"]
        set refBa [atomselect top "same residue as index $ba"]
        set refAcIndex [lindex [$refAc list] 0]
        set refBaIndex [lindex [$refBa list] 0]


        if [info exists finalPairs($refAcIndex,$refBaIndex)] {
          unset potPairs($ac,$ba)
        } else {
          set acName [lindex [$refAc get resname] 0]
          set acId [lindex [$refAc get resid] 0]
          set acSegname [lindex [$refAc get segname] 0]
          set baName [lindex [$refBa get resname] 0]
          set baId [lindex [$refBa get resid] 0]
          set baSegname [lindex [$refBa get segname] 0]
          if {$acSegname ne $baSegname} {
            set finalPairs($refAcIndex,$refBaIndex) 1
            set out_line [join "PRINTING: $acSegname $acName $acId - $baSegname $baName $baId"]
            puts "$out_line"
          }
        }
        $refAc delete
        $refBa delete
    }
    set res1 [llength [array names finalPairs]]

    array unset potPairs
    array unset finalPairs
    
    set sel1 [atomselect top "(protein and basic and nitrogen and not backbone) and segname D$seg and resid $sel_input1" frame $nn]
    set sel2 [atomselect top "(protein and acidic and oxygen and not backbone) and not segname D$seg and resid $sel_input2" frame $nn]
    set tmpList [measure contacts $ondist $sel1 $sel2]
    $sel1 delete
    $sel2 delete
    
    foreach i [lindex $tmpList 0] j [lindex $tmpList 1] {
        set potPairs($i,$j) 1
    }

    foreach pair [array names potPairs] {
    
        foreach { ac ba } [split $pair ,] break

        set refAc [atomselect top "same residue as index $ac"]
        set refBa [atomselect top "same residue as index $ba"]
        set refAcIndex [lindex [$refAc list] 0]
        set refBaIndex [lindex [$refBa list] 0]


        if [info exists finalPairs($refAcIndex,$refBaIndex)] {
          unset potPairs($ac,$ba)
        } else {
          set acName [lindex [$refAc get resname] 0]
          set acId [lindex [$refAc get resid] 0]
          set acSegname [lindex [$refAc get segname] 0]
          set baName [lindex [$refBa get resname] 0]
          set baId [lindex [$refBa get resid] 0]
          set baSegname [lindex [$refBa get segname] 0]
          if {$acSegname ne $baSegname} {
            set finalPairs($refAcIndex,$refBaIndex) 1
            set out_line [join "PRINTING: $acSegname $acName $acId - $baSegname $baName $baId"]
            puts "$out_line"
          }
        }
        $refAc delete
        $refBa delete
    }
    set res2 [llength [array names finalPairs]]

    return [format "%d" [expr $res1 + $res2]]
}