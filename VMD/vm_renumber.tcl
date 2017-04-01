mol new ../combine_dcd/initial.psf
mol addfile ../combine_dcd/initial.pdb
set sel1 [atomselect top all]
$sel1 set segname L1
$sel1 set chain L
set oldNumbering [$sel1 get residue]
set newNumbering {}
foreach resid $oldNumbering {
    lappend newNumbering [expr {$resid + 1}]
}
$sel1 set resid $newNumbering

$sel1 writepdb renumbered.pdb
$sel1 writepsf renumbered.psf