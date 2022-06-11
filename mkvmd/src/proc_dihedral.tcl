# Function as it is named
# Only for the function dihedral
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
  return [expr $sign * 57.2958 * acos($dotprod / ($amag * $bmag))]
}

# Function as it is named
# dihedral <pos_a1> <pos_a2> <pos_a3> <pos_a4>
# a1
#   \
#    (a2,a3)
#           \
#            a4
proc dihedral { a1 a2 a3 a4 } {
  set r1 [vecsub $a1 $a2]
  set r2 [vecsub $a3 $a2]
  set r3 [vecscale $r2 -1]
  set r4 [vecsub $a4 $a3]

  set n1 [veccross $r1 $r2]
  set n2 [veccross $r3 $r4]
  
  return [signed_angle $n1 $n2 $r2]
}