set psfname [lindex $argv 0]
set dcdname [lindex $argv 1]

set molnum [mol new $psfname waitfor all]
mol addfile $dcdname waitfor all

package require pbctools
#pbc set {121 121 356}
pbc unwrap -all

animate write dcd $dcdname-unwrap.dcd $molnum

quit
