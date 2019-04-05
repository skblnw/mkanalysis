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
    # Therefore you could define the box size manually according to the newest xsc
    set vec_def {100 100 150}
}

# In fact, you DO NOT have to define pbc box if you are only wrapping dcd. Box size already exists.
# Wrapping the trajectory according to the pbc size
set PREFIX noW
set PSFNAME initial-$PREFIX.psf
#set PDBNAME initial-$PREFIX.pdb

set PREFIX_DCD npt-noW-pf1000ps-1
set DCDNAME ${PREFIX_DCD}.dcd

set mol [mol new $PSFNAME type psf waitfor all]
#mol addfile $PDBNAME type pdb waitfor all
mol addfile $DCDNAME type dcd waitfor all

package require pbctools
# pbc set $vec_init
# pbc set $vec_def
# Sometimes it is more straight-forward to read a xsc file for pbc size
# pbc read
# Sometimes you have to unwrap the trajectory first if the protein jumped too much
# pbc unwrap -all
# In fact, you DO NOT have to define pbc box if you are only wrapping dcd. Box size already exists.
pbc wrap -first 1 -last last -center com -centersel "protein" -compound residue

set comseltext "protein"

set frames [molinfo $mol get numframes]
set all [atomselect $mol all]
set ref [atomselect $mol ($comseltext) frame 0]
set comsel [atomselect $mol ($comseltext)]

# Align to first frame
if {1} {
puts "Info) Done wrapping, now align selection to the first frame"
puts "Info) You have to figure out by yourself whether the selection has been moved to the origin. Default: NO"
for {set i 0} {$i < $frames} {incr i} {
  $all frame $i
  $all update
  $comsel frame $i
  $comsel update
  $all move [measure fit $comsel $ref]

  # Move to center
  # You may decide here whether to move the protein to the origin or not
  if {0} {
      set com [measure center $comsel weight mass ]
      # set lat [molinfo $mol get {a b c}]
      # set shift [vecsub [vecscale $lat 0.5] $com]
      # $all moveby $shift
      $all moveby [vecinvert $com]
  }
}
}

puts "Info) Done, now write wrapped dcd"

animate write dcd ${PREFIX_DCD}-wrapped.dcd $mol
quit
