proc measureSasa { nn probe_size resid } {
    set rest [atomselect top "segname PROA and resid $resid and not name C N O CA HA HN" frame $nn]
    set selprot [atomselect top "segname PROA" frame $nn]
    set sasa [measure sasa $probe_size $selprot -restrict $rest]
    if { $sasa != 0 } {
        return 1
    } else {
        return 0
    }
}

# /------------------/
# /     Main Body    /
# /------------------/

# !!!Important!!!
# Deleting existing files as we APPEND instead of trashing and opening new files
# Make sure you will delete all the existing files
# !!!Important!!!
# eval file delete [glob output/*.dat]
# set OUTPUT_DIR [exec date +%Y%m%d%H%M%S]
# exec mkdir -p $OUTPUT_DIR

# mol new $PDB waitfor all
# mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]


set selprot [atomselect top "segname PROA"]
puts "mkvmd> These residues will be scanned: "

set sel [atomselect top "{segname PROA and not name C N O CA HA HN} and within 4 of segname ANTI "]
set tmplist []
foreach resid [lsort -unique -integer [$sel get resid]] { 
    set rest [atomselect top "segname PROA and resid $resid"]
    set sasa [measure sasa 5.0 $selprot -restrict $rest]
    if { $sasa == 0 } { 
        lappend tmplist $resid 
    }
}
set reslist []
foreach resid $tmplist {  
    set rest [atomselect top "segname PROA and resid $resid"]
    set sasa [measure sasa 1.4 $selprot -restrict $rest]
    if { $sasa != 0 } { 
        puts -nonewline "$resid "
        lappend reslist $resid 
    }
}
puts ""
puts "mkvmd> [llength $reslist] residues in total"

puts "mkvmd> Computing something..."
for {set nn 0} {$nn < $total_frame} {incr nn} {
    set nframe [expr $nn + 0]
    puts "Frame: $nframe"
    $selprot frame $nn
    $selprot set user 10

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
    # This determines <number of output files>
    foreach {sel_input1} {1} {sel_input2} {1} {
      # Output file name
      # set outf [open ${OUTPUT} "a"]
      # Write TIME at the very first of a line
      # set out_line [format "%d" \$nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
      # This determines <number of columns in one output file>

      # if { \$nn == 0 } {
      #   foreach resid \$reslist {
      #       lappend out_line \$resid
      #   }
      #   puts \$outf "\$out_line"
      #   set out_line [format "%d" \$nframe]
      # }

      foreach resid $reslist {
        # Call calc funtion you like
        # set out_line [format "%.1f" \$probe_size]
        set count 0
        foreach probe_size {2.5 3.0 3.5 4.0} {
            incr count [measureSasa $nn $probe_size $resid]
        }
        # lappend out_line \$count
        puts -nonewline "$count "
        set sel [atomselect top "segname PROA and resid $resid" frame $nn]
        $sel set user $count
      }
      # puts \$outf "\$out_line"

      # Write to file
      # Remember to close the file
      # close \$outf
    }
}

puts ""
