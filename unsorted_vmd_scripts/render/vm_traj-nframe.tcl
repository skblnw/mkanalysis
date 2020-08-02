proc enabletrace {} { 
  global vmd_frame; 
  trace variable vmd_frame([molinfo top]) w drawcounter 
} 
proc disabletrace {} { 
  global vmd_frame; 
  trace vdelete vmd_frame([molinfo top]) w drawcounter 
} 
proc drawcounter { name element op } { 
  global vmd_frame; 
  draw delete all 
  # puts "callback!" 
  draw color black 
  set nsperframe .1
  set time [format "%3.2f ns" [expr $vmd_frame([molinfo top]) * $nsperframe]] 
  draw text {0 0 0} "$time" size 2 thickness 5 
}