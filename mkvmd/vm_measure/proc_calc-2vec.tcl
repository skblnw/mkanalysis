source proc_angle.tcl

# Calc function of angles between 2 principle axes
proc calc_2vec {nn vec1 vec2} {
        set ang [angle $vec1 $vec2]

        set res  [format "%.2f" $ang]
        return $res
}
