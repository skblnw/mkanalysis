proc findOverlap { nn ref dist tt} {
    set sel [atomselect top "name CA and same residue as {chain B and within $dist of chain A}" frame $nn]

    foreach i [$sel get resid] {
        set tarRes($i) 1
    }

    set cc 0
    foreach res [$ref get resid] {
        if [info exist tarRes($res)] {
            incr cc
        }
    }
    
    return [format "%4.2f" [expr {$cc / ($tt + 0.0)}]]
}