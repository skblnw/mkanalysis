mol new ../../../ionized.pdb
set sel [atomselect top protein]
$sel writepdb first.pdb
mol delete all

package require autopsf
mol new first.pdb
autopsf -mol top -top top_all36_prot.rtf
