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

# Calc function of angles between 2 principle axes
proc calc_2pa {nn Ip1 Ip2} {
        set ang [angle $Ip1 $Ip2]
        set res  [format "%.2f" $ang]

        return $res
}

# Load packages for calculating principle axis (if needed)
package require Orient 
namespace import Orient::orient

# set sel [atomselect top "protein and segname P2 and name CA"]
foreach ii {1} {
	set sel [atomselect top "protein and segname P$ii and name CA"]
	set I [draw principalaxes $sel]
	set sel [atomselect top "protein and segname P$ii and resid 345 to 361 and name CA"]
	set I [draw principalaxes $sel]
}