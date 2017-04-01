puts "Info) Loading ionized.pdb for wrapping box size"
mol new ../../../ionized.psf waitfor all
mol addfile ../../../ionized.pdb waitfor all

set sel [atomselect top {segname "WT.*"}]
set minmax [measure minmax $sel]
set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
$sel delete
mol delete all
puts "Info) Box size got, delete top molecule, loading trajectory"

set DCDNAME npt-pf1000ps-307.dcd

set mol [mol new initial-noW.psf waitfor all]
mol addfile $DCDNAME waitfor all

package require pbctools
pbc set $vec
pbc unwrap -all
pbc wrap -all -center com -centersel "protein"
puts "Info) Done wrapping, now move selection to center"

set comseltext "protein"

set frames [molinfo $mol get numframes]
set all [atomselect $mol all]
set comsel [atomselect $mol ($comseltext)]

for {set i 0} {$i < $frames} {incr i} {
  $all frame $i
  $all update
  $comsel frame $i
  $comsel update

  set com [measure center $comsel weight mass ]
  set lat [molinfo $mol get {a b c}]
  set shift [vecsub [vecscale $lat 0.5] $com]

  $all moveby $shift
}
puts "Info) Done, now write wrapped dcd"

animate write dcd wrapped.dcd $mol
