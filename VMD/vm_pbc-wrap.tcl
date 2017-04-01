if {0} {
    # Obtaining pbc size from water molecules in ionized.pdb
    # Note that the box size could be very different after equilibrations
    puts "Info) Loading ionized.pdb for wrapping box size"
    mol new ../../../ionized.psf waitfor all
    mol addfile ../../../ionized.pdb waitfor all
    
    set sel [atomselect top {segname "WT.*"}]
    set minmax [measure minmax $sel]
    set vec_init [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    $sel delete
    mol delete all
    puts "Info) Box size got, delete top molecule, loading trajectory"
}

if {0} {
    # Therefore you could (and are recommended to) define the box size manually according to the newest xsc
    set vec_def {100 100 150}
}

# Wrapping the trajectory according to the pbc size
set PSFNAME initial-noW.psf
set PDBNAME initial-noW.pdb
set DCDNAME Bigdcd/npt-pf200ps-3-full.dcd

set mol [mol new $PSFNAME type psf waitfor all]
mol addfile $PDBNAME type pdb waitfor all
mol addfile $DCDNAME type dcd waitfor all

package require pbctools
# pbc set $vec_init
# pbc set $vec_def
# Sometimes it is more straight-forward to read a xsc file for pbc size
# pbc read
# Sometimes you have to unwrap the trajectory first if the protein jumped too much
# pbc unwrap -all
# In fact, you do NOT have to define pbc box if you are only wrapping dcd. Box size already exists.
pbc wrap -first 1 -last last -center com -centersel "segname P1" -compound fragment

set comseltext "protein"

set frames [molinfo $mol get numframes]
set all [atomselect $mol all]
set ref [atomselect $mol ($comseltext) frame 0]
set comsel [atomselect $mol ($comseltext)]

# Align to first frame
if {1} {
puts "Info) Done wrapping, now align selection to the first frame"
for {set i 1} {$i < $frames} {incr i} {
  $all frame $i
  $all update
  $comsel frame $i
  $comsel update

  $all move [measure fit $comsel $ref]
}
}

# Move to center
if {0} {
puts "Info) Done wrapping, now move selection to center"
for {set i 1} {$i < $frames} {incr i} {
  $all frame $i
  $all update
  $comsel frame $i
  $comsel update

  set com [measure center $comsel weight mass ]
  set lat [molinfo $mol get {a b c}]
  set shift [vecsub [vecscale $lat 0.5] $com]

  $all moveby $shift
}
}
puts "Info) Done, now write wrapped dcd"

animate write dcd Bigdcd/wrapped-pf200ps.dcd $mol
quit
