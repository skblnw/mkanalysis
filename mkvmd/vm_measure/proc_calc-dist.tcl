# Calc function as it is named
proc calc_dist {nn sel_input1 sel_input2} {
    set sel1 [atomselect top "protein and resid ${sel_input1} and name CA" frame $nn]
    set sel2 [atomselect top "not resname "D.*" and resid ${sel_input2} and name CA" frame $nn]
    
    set dist [vecdist [measure center $sel1 weight mass] [measure center $sel2 weight mass]]
    set res [format "%.2f" $dist]

    $sel1 delete
    $sel2 delete

    return $res
}