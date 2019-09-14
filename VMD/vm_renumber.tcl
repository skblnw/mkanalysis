foreach ii {A B C D} {
    mol new $ii.pdb
    set sel [atomselect top all]
    set oldNumbering [$sel get residue]
    set newNumbering {}
    foreach resid $oldNumbering {
        lappend newNumbering [expr {$resid + 1}]
    }
    $sel set resid $newNumbering

    $sel writepdb $ii.pdb
}

