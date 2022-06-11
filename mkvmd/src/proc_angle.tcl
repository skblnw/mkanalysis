# Function as it is named
# This is the 'real' function for calculating angles (unless I find the signed_angle at anytime usedful)
# angle <vec1> <vec2>
#      o   o
#  vec1 \ / vec2
#        o
proc angle { a b } {
  # get Pi 
  global M_PI
   # Angle between two vectors 
  set cosine [expr [vecdot $a $b] / ( [veclength $a] * [veclength $b])]
  return [expr acos($cosine)*(180.0/$M_PI)]
}