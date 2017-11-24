# Find contacts within "dist"
# Add atom indexes to selections to make the script compute results for atom pairs
proc findContacts { nn dist } {
    set sel1 [atomselect top "protein and segname P1 and resid 250 to 361 and name CA" frame $nn]
    foreach id [$sel1 get index] {
        set sel [atomselect top "index $id" frame $nn]
        set selres [atomselect top "same residue as index $id" frame $nn]
        $sel set {x y z} [list [measure center $selres weight mass]]
    }
    set sel2 [atomselect top "resname PIP2 and resid 95 100 and name P" frame $nn]
    set tmpList [measure contacts $dist $sel1 $sel2]
    $sel1 delete
    $sel2 delete
    
    foreach i [lindex $tmpList 0] j [lindex $tmpList 1] {
        set potPairs($i,$j) 1
    }

    foreach pair [array names potPairs] {
        foreach { ac ba } [split $pair ,] break

        set refAc [atomselect top "index $ac"]
        set refBa [atomselect top "index $ba"]
        set acName [lindex [$refAc get resname] 0]
        set acId [lindex [$refAc get resid] 0]
        set acSegname [lindex [$refAc get segname] 0]
        set baName [lindex [$refBa get resname] 0]
        set baId [lindex [$refBa get resid] 0]
        set baSegname [lindex [$refBa get segname] 0]
        # set finalPairs($refAcIndex,$refBaIndex) 1
        set out_line [join "FRAME $nn: $acSegname $acName $acId - $baSegname $baName $baId"]
        puts "$out_line"
        $refAc delete
        $refBa delete
    }
    set res [format "%d" [llength [array names potPairs]]]
    return $res
}
