set filename ""

package require Orient
namespace import Orient::orient

set sel [atomselect top "protein and segname P2 and resid 250 to 361 and name CA"]
set selall [atomselect top "protein"]

draw delete all
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$selall move $A