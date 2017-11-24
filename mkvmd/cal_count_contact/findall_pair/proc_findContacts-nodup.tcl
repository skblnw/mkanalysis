# Find contacts within "dist"
# Add atom indexes to selections to make the script compute results for atom pairs
proc findContacts_nodup { nn seg dist } {
    set sel1 [atomselect top "protein and segname D$seg and name CA" frame $nn]
    set sel2 [atomselect top "protein and not segname D$seg and name CA" frame $nn]
    set tmpList [measure contacts $dist $sel1 $sel2]
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
            set out_line [join "FOUND: $acSegname $acName $acId - $baSegname $baName $baId"]
            puts "$out_line"
          }
        }
        $refAc delete
        $refBa delete
    }
    set res [format "%d" [llength [array names finalPairs]]]
    return $res
}