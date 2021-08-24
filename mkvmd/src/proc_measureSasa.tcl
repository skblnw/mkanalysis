proc calSasa { nn } {
    set sel1 [atomselect top "chain A" frame $nn]
    set sel2 [atomselect top "chain B" frame $nn]
    set selall [atomselect top "chain A B" frame $nn]
    
    return [format "%.2f" [expr ([measure sasa 1.4 $selall] - [measure sasa 1.4 $sel1] - [measure sasa 1.4 $sel2]) / 2.0]]
}