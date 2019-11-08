## Transpose a 2D table (list of lists). From [[http://wiki.tcl.tk/2748]]
proc transpose {matrix} {
    set cmd list
    set i -1
    foreach col [lindex $matrix 0] {append cmd " \$[incr i]"}
    foreach row $matrix {
        set i -1
        foreach col $row {lappend [incr i] $col}
    }
    eval $cmd
}

## Prepare the "reference" list of native contacts, e.g. from the
# crystal structure. Return value: a list of native contact pairs
# (only useful to be passed as an argument to measureNativeContacts,
# or to get its length).
#
# *Note.* The function is subject to change. Later can save sizes of
# atomselections (for safety checking) and/or cutoff. 
proc prepareNativeContacts { cutoff sel1 {sel2 0} } {
    if { $sel2 == 0 } { set sel2 $sel1 }
    return [ transpose [ measure contacts $cutoff $sel1 $sel2 ] ]
}

proc measureNativeContacts { nclist cutoff sel1 {sel2 0} } {
    set n 0
    if { $sel2 == 0 } { set sel2 $sel1 }
    # current contacts
    set cc [ transpose [ measure contacts $cutoff $sel1 $sel2 ] ]
    # check if current is in native set (intersect)
    # iterate over current pairs cp
    foreach cp $cc {
        if {    [lsearch $nclist $cp] != -1 ||
            [lsearch $nclist [lreverse $cp]] != -1 } {incr n }
    }
    return $n
}
