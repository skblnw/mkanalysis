# Function as it is named
# Solely for the function dihedral
# signed_angle <pos_a> <pos_b> <pos_c>
# a   c
#  \ /
#   b
proc signed_angle { a b c } {
  set amag [veclength $a]
  set bmag [veclength $b]
  set dotprod [vecdot $a $b]
  set crossp [veccross $a $b]
  set sign [vecdot $crossp $c]
  if { $sign < 0 } {
    set sign -1
  } else {
    set sign 1
  }
  return [expr $sign * (180.0/$M_PI) * acos($dotprod / ($amag * $bmag))]
}
