set psfname [lindex $argv 0]
set dcdname [lindex $argv 1]
set mapname [lindex $argv 2]
set mapres [lindex $argv 3]
set outname [lindex $argv 4]
set thres [lindex $argv 5]

set molnum [mol new $psfname waitfor all]
mol addfile $dcdname waitfor all

package require mdff
mdff check -mol 0 -frames all -ccc -map $mapname -res $mapres -cccfile $outname_ccc_global.dat
mdff check -mol 0 -frames all -ccc -map $mapname -res $mapres -cccfile $outname_ccc_local.dat -threshold $thres

quit
