# Function as it is named
# This function returns the projection of vector r onto plane defined by vector a and b

proc proj2plane { r a b } {
    set n [veccross $a $b]
    return [vecsub $r [vecscale $n [vecdot $r $n]]]
}