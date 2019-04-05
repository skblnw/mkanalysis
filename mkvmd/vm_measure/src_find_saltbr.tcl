package provide find_saltbr 1.1

namespace eval ::find_saltbr:: {
  namespace export find_saltbr

  variable defaultONDist 3.2
  variable defaultFrames "all"
  variable defaultOutdir
  variable defaultLogFile ""
  variable defaultUpdateSel 1
  variable debug 0
  variable currentMol none
  variable atomselectText "protein"
  variable statusMsg ""
}

proc ::saltbr::saltbr_usage { } {

  variable defaultCOMDist
  variable defaultONDist
  variable defaultWriteAll
  variable defaultFrames

  puts "Usage: saltbr -sel <atom selection> <option1> <option2> ..."
  #puts "  -sel <atom selection> (default: \[atomselect top protein\])"
  puts "Options:"
  puts "  -upsel <yes|no> (update atom selections every frame? default: yes)"
  puts "  -frames <begin:end> or <begin:step:end> or all or now (default: $defaultFrames)"
  puts "  -ondist <cutoff distance between oxygen and nitrogen atoms> (default: $defaultONDist)"
  puts "  -outdir <output directory> (default: current)"
  puts "  -log <log filename> (default: none)"
  return
}

proc ::saltbr::saltbr { args } {

  variable defaultONDist
  variable defaultFrames
  variable defaultFrames
  variable defaultUpdateSel
  variable currentMol
  variable atomselectText
  variable debug
  variable log
  variable statusMsg
  
  variable defaultOutdir [pwd]

  set nargs [llength $args]
  if { $nargs == 0 || $nargs % 2 } {
    if { $nargs == 0 } {
      saltbr_usage
      error ""
    }
    if { $nargs % 2 } {
      saltbr_usage
        error "error: odd number of arguments $args"
    }
  }

  foreach {name val} $args {
    switch -- $name {
      -sel { set arg(sel) $val }
      -upsel { set arg(upsel) $val }
      -frames { set arg(frames) $val }
      -ondist { set arg(ondist) $val }
      -outdir { set arg(outdir) $val }
      -log { set arg(log) $val }
      -debug { set arg(debug) $val }
      default { error "unkown argument: $name $val" }
    }
  }

  # debug flag
  if [info exists arg(debug)] {
      set debug 1
  }

  # outdir
  if [info exists arg(outdir)] {
    set outdir $arg(outdir)
  } else {
    set outdir $defaultOutdir
  }
  if { ![file isdirectory $outdir] } {
    error "$outdir is not a directory."
  }

  # log file
  if { [info exists arg(log)] && $arg(log) != "" } {
    set log [open [file join $outdir $arg(log)] w]
  } else {
    set log "stdout"
  }

  # get selection
  if [info exists arg(sel)] {
    set sel $arg(sel)
    set molid [$sel molid]
  } elseif $gui {
    if { $currentMol == "none" } {
      error "No molecules were found."
    } else {
      set molid $currentMol
      set sel [atomselect $currentMol $atomselectText]
    }
  } else {
    saltbr_usage
    error "No atomselection was given."
  }

  # update selections?
  if [info exists arg(upsel)] {
    if { $arg(upsel) == "no" || $arg(upsel) == 0 } {
      set updateSel 0
    } elseif { $arg(upsel) == "yes" || $arg(upsel) == 1 } {
      set updateSel 1
    } else {
      error "error: bad argument for option -upsel $arg(upsel): acceptable arguments are 'yes' or 'no'"
    }
  } else {
    set updateSel $defaultUpdateSel
  }

  # get frames
  set nowframe [molinfo $molid get frame]
  set lastframe [expr [molinfo $molid get numframes] - 1]
  if { ! [info exists arg(frames)] } { set arg(frames) $defaultFrames }
  if [info exists arg(frames)] {
    set fl [split $arg(frames) :]
    switch -- [llength $fl] {
      1 {
        switch -- $fl {
          all {
            set frames_begin 0
            set frames_end $lastframe
          }
          now {
            set frames_begin $nowframe
          }
          last {
            set frames_begin $lastframe
          }
          default {
            set frames_begin $fl
          }
        }
      }
      2 {
        set frames_begin [lindex $fl 0]
        set frames_end [lindex $fl 1]
      }
      3 {
        set frames_begin [lindex $fl 0]
        set frames_step [lindex $fl 1]
        set frames_end [lindex $fl 2]
      }
      default { error "bad -frames arg: $arg(frames)" }
    }
  } else {
    set frames_begin 0
  }
  if { ! [info exists frames_step] } { set frames_step 1 }
  if { ! [info exists frames_end] } { set frames_end $lastframe }
    switch -- $frames_end {
      end - last { set frames_end $lastframe }
  }
  if { [ catch {
    if { $frames_begin < 0 } {
      set frames_begin [expr $lastframe + 1 + $frames_begin]
    }
    if { $frames_end < 0 } {
      set frames_end [expr $lastframe + 1 + $frames_end]
    }
    if { ! ( [string is integer $frames_begin] && \
      ( $frames_begin >= 0 ) && ( $frames_begin <= $lastframe ) && \
    [string is integer $frames_end] && \
      ( $frames_end >= 0 ) && ( $frames_end <= $lastframe ) && \
      ( $frames_begin <= $frames_end ) && \
      [string is integer $frames_step] && ( $frames_step > 0 ) ) } {
        error
      }
  } ok ] } { error "bad -frames arg: $arg(frames)" }
  if $debug {
    puts $log "frames_begin: $frames_begin"
    puts $log "frames_step: $frames_step"
    puts $log "frames_end: $frames_end"
    flush $log
  }

  # get ONDIst
  if [info exists arg(ondist)] {
    set ONDist $arg(ondist)
  } else {
    set ONDist $defaultONDist
  }

  # print name, version and date of plugin
  puts $log "Salt Bridges Plugin, Version 1.0"
  puts $log "[clock format [clock scan now]]\n"
  puts $log "Parameters used in the calculation of salt bridges:"
  puts $log "- Atomselection: [$sel text]"
  if $updateSel {
    puts $log "- Update selections every frame: yes"
  } else {
    puts $log "- Update selections every frame: no"
  }
  puts $log "- Initial frame: $frames_begin"
  puts $log "- Frame step: $frames_step"
  puts $log "- Final frame: $frames_end"
  puts $log "- Oxygen-nitrogen cut-off: $ONDist"
  puts $log "- Center of mass cut-off: $COMDist"
  if $writefiles {
    puts $log "- Write a file for each salt bridge: yes"
  } else {
    puts $log "- Write a file for each salt bridge: no"
  }
  puts $log ""
  flush $log
  

  for { set nn $frames_begin } { $nn <= $frames_end } { incr nn 1 } {
    foreach seltext {P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P15 P16 P17 P18 P19 P20 P21 P22 P23 P24}
      # pairs is an associative array containing the salt bridges
      set numpairs [findSaltBridges2 $nn $seltext $ONDist]
    }
  }

  # delete the selection if it was created here
  if { ![info exists arg(sel)] } {
    $sel delete
  }

  if { $log != "stdout" } {
      close $log
  }

  set statusMsg "Done."
  update

  return $numpairs

}

proc ::saltbr::findSaltBridges { nn seltext ondist} {
    set sel1 [atomselect top "(protein and acidic and oxygen and not backbone) and $seltext" frame $nn]
    set sel2 [atomselect top "(protein and basic and nitrogen and not backbone) and not $seltext" frame $nn]
    set tmpList [measure contacts $ondist $sel1 $sel2]
    $sel1 delete
    $sel2 delete
    
    foreach i [lindex $tmpList 0] j [lindex $tmpList 1] {
        set potPairs($i,$j) 1
    }

    foreach pair [array names potPairs] {
    
        foreach { ac ba } [split $pair ,] break

        set refAc [atomselect top "same residue as index $ac"]
        set refBa [atomselect top "same residue as index $ba"]
        set refAcIndex [lindex [$refAc list] 0]
        set refBaIndex [lindex [$refBa list] 0]


        if [info exists finalPairs($refAcIndex,$refBaIndex)] {
          unset potPairs($ac,$ba)
        } else {
          set acName [lindex [$refAc get resname] 0]
          set acId [lindex [$refAc get resid] 0]
          set acSegname [lindex [$refAc get segname] 0]
          set baName [lindex [$refBa get resname] 0]
          set baId [lindex [$refBa get resid] 0]
          set baSegname [lindex [$refBa get segname] 0]
          if {$acSegname ne $baSegname} {
            set finalPairs($refAcIndex,$refBaIndex) 1
            set out_line [join "$acName$acId-segname$acSegname $baName$baId-segname$baSegname"]
            puts "$out_line"
          }
        }
        $refAc delete
        $refBa delete
    }
  }

    set res [format "%d" [llength [array names finalPairs]]]
    return $res
}

proc ::saltbr::findSaltBridges2 { selection updateSel COMDist ONDist frames_begin frames_step frames_end pairsName idpairsName } {

  variable debug
  variable statusMsg
  variable log

  if $debug {
    puts $log "updateSel = $updateSel"
  }

  upvar $pairsName finalPairs
  upvar $idpairsName idPairs

  set molid [$selection molid]
  set seltext [$selection text]

  set acsel [atomselect $molid "(protein and acidic and oxygen and not backbone) and $seltext"]
  set basel [atomselect $molid "(protein and basic and nitrogen and not backbone) and $seltext"]

  if { [$acsel num] == 0 || [$basel num] == 0 } {
    if { [$acsel num] == 0 } {
      set errMsg "No oxygens of acidic amino acid residues were found in the given selection."
    }
    if { [$basel num] == 0 } {
      append errMsg "\nNo nitrogens of basic amino acid residues were found in the given selection."
    }
    if { $log != "stdout" } {
      close $log
    }
    error $errMsg
  }

  set statusMsg "Searching for ion-pairs with oxygen-nitrogen distance\nwithin $ONDist Angstroms in the selected frames... "
  set statusMsg2 "Searching for ion-pairs with oxygen-nitrogen distance within $ONDist Angstroms in the selected frames... "
  update
  puts -nonewline $log $statusMsg2
  flush $log

  for { set f $frames_begin } { $f <= $frames_end } { incr f $frames_step } {
    $acsel frame $f
    $basel frame $f
    if $updateSel {
      $acsel update
      $basel update
    }
    set tmpList [measure contacts $ONDist $acsel $basel] 
      foreach i [lindex $tmpList 0] j [lindex $tmpList 1] {
        set potPairs($i,$j) 1
      }
  }
  append statusMsg "Done."
  puts $log "Done."
  flush $log

  # Remove redundancies in the list of salt bridges
  set statusMsg "Removing redundancies in the ion-pairs found... "
  update
  puts -nonewline $log $statusMsg
  flush $log
  foreach pair [array names potPairs] {
    
    foreach { ac ba } [split $pair ,] break

    set refAc [atomselect $molid "same residue as index $ac"]
    set refBa [atomselect $molid "same residue as index $ba"]
    set refAcIndex [lindex [$refAc list] 0]
    set refBaIndex [lindex [$refBa list] 0]
    $refAc delete
    $refBa delete

    if [info exists refPotPairs($refAcIndex,$refBaIndex)] {
      unset potPairs($ac,$ba)
    } else {
      set refPotPairs($refAcIndex,$refBaIndex) 1
    }

  }
  append statusMsg "Done."
  puts $log "Done."
  flush $log

  $acsel delete
  $basel delete


  if { $COMDist == "none" } {
    foreach key [array names potPairs] {
      set finalPairs($key) $potPairs($key)
    }
  } else {
    set statusMsg "Selecting ion-pairs whose side chains' centers of mass\nare within $COMDist Angstroms... "
    set statusMsg2 "Selecting ion-pairs whose side chains' centers of mass are within $COMDist Angstroms... "
    update
    puts -nonewline $log $statusMsg2
    flush $log
    
    set aclist [list]
    set balist [list]
    foreach pair [array names potPairs] {
  
      foreach {ac ba} [split $pair ,] break
      
      # select heavy atoms of the side chains
      set acsel [atomselect $molid "not backbone and noh and same residue as index $ac"]
      set basel [atomselect $molid "not backbone and noh and same residue as index $ba"]
            
      for { set f $frames_begin } { $f <= $frames_end } { incr f $frames_step } {
        $acsel frame $f
        $basel frame $f
        if $updateSel {
          $acsel update
          $basel update
        }
    
        # find the center of mass of each side-chain
        set accenter [measure center $acsel weight mass]
        set bacenter [measure center $basel weight mass]

        if { [veclength [vecsub $accenter $bacenter]] <= $COMDist } {
    lappend aclist $ac
          lappend balist $ba
          break
        }
      }
            
      $acsel delete
      $basel delete
    }

    foreach ac $aclist ba $balist {
      set finalPairs($ac,$ba) 1
    }

    append statusMsg "Done."
    puts $log "Done."
    flush $log

  }

  # extract identification of each pairs and store in idPairs
  set statusMsg "Extracting identification of each salt bridge... "
  update
  puts -nonewline $log $statusMsg
  flush $log

  set useChain 0
  set useSegname 0
  set extractId 1

  while { $extractId != 0 } {

    set extractId 0

    foreach pair [array names finalPairs] {
  
      foreach {ac ba} [split $pair ,] break
  
      # select heavy atoms of the side chains
      set acsel [atomselect $molid "not backbone and noh and same residue as index $ac"]
      set basel [atomselect $molid "not backbone and noh and same residue as index $ba"]
      set acseloxy [atomselect $molid "not backbone and oxygen and same residue as index $ac"]
      set baselnit [atomselect $molid "not backbone and nitrogen and same residue as index $ba"]
          
      # get resid, chain, and segname
      set acresid [lsort -unique [$acsel get resid]]
      set acresname [lsort -unique [$acsel get resname]]
      set acchain [lsort -unique [$acsel get chain]]
      set acsegname [lsort -unique [$acsel get segname]]
      set baresid [lsort -unique [$basel get resid]]
      set baresname [lsort -unique [$basel get resname]]
      set bachain [lsort -unique [$basel get chain]]
      set basegname [lsort -unique [$basel get segname]]
         
      set seltext "protein and not backbone and noh"
      set acid "$acresname$acresid"
      set baid "$baresname$baresid"

      # use chain in the identification?
      set acseltest [atomselect $molid "$seltext and resid $acresid"]
      if { $useChain != 0 && $acchain != 0 } {
        append acid "_chain$acchain"
      } elseif { [$acsel num] != [$acseltest num] && $acchain != 0 } {
        set useChain 1
        set extractId 1
        break
      }
      $acseltest delete

      # use segname in the identification?
      set acseltest [atomselect $molid "$seltext and resid $acresid and chain $acchain"]
      if { $useSegname != 0 && $acsegname != 0 && $acsegname != "{}" } {
        append acid "_segname$acsegname"
      } elseif { [$acsel num] != [$acseltest num] && $acsegname != 0 && $acsegname != "{}" } {
        set useSegname 1
        set extractId 1
        break
      }
      $acseltest delete

      # is this enough to identify the residue?
      if { $useChain != 0 && $useSegname != 0 } {
        set acseltest [atomselect $molid "$seltext and resid $acresid and chain $acchain and segname $acsegname"]
      } elseif { $useChain != 0 } {
        set acseltest [atomselect $molid "$seltext and resid $acresid and chain $acchain"]
      } elseif { $useSegname != 0 } {
        set acseltest [atomselect $molid "$seltext and resid $acresid and segname $acsegname"]
      } else {
        set acseltest [atomselect $molid "$seltext and resid $acresid"]
      }
      if { [$acsel num] != [$acseltest num] } {
        set statusMsg "\nWarning: the identification $acid is not unique."
        update
        puts $log $statusMsg
        flush $log
      }
      $acseltest delete

      # use chain in the identification?
      set baseltest [atomselect $molid "$seltext and resid $baresid"]
      if { $useChain != 0 && $bachain != 0 } {
        append baid "_chain$bachain"
      } elseif { [$basel num] != [$baseltest num] && $bachain != 0 } {
        set useChain 1
        set extractId 1
        break
      }
      $baseltest delete

      # use segname in the identification?
      set baseltest [atomselect $molid "$seltext and resid $baresid and chain $bachain"]
      if { $useSegname != 0 && $basegname != 0 && $basegname != "{}" } {
        append baid "_segname$basegname"
      } elseif { [$basel num] != [$baseltest num] && $basegname != 0 && $basegname != "{}" } {
        set useSegname 1
        set extractId 1
        break
      }
      $baseltest delete

      # is this enough to identify the residue?
      if { $useChain != 0 && $useSegname != 0 } {
        set baseltest [atomselect $molid "$seltext and resid $baresid and chain $bachain and segname $basegname"]
      } elseif { $useChain != 0 } {
        set baseltest [atomselect $molid "$seltext and resid $baresid and chain $bachain"]
      } elseif { $useSegname != 0 } {
        set baseltest [atomselect $molid "$seltext and resid $baresid and segname $basegname"]
      } else {
        set baseltest [atomselect $molid "$seltext and resid $baresid"]
      }
      if { [$basel num] != [$baseltest num] } {
        set statusMsg "\nWarning: the identification $baid is not unique."
        update
        puts $log $statusMsg
        flush $log
      }
      $baseltest delete
  
      set idPairs($pair) "$acid-$baid"
  
    }

  }

  append statusMsg "Done."
  puts $log "Done."
  flush $log

  set numpairs [llength [array names finalPairs]]
  set statusMsg "Found $numpairs salt bridges."
  update
  puts $log "$statusMsg\n"
  flush $log

  foreach pair [array names idPairs] {
      puts $log $idPairs($pair)
  }
  puts $log ""
  flush $log

  return $numpairs

}