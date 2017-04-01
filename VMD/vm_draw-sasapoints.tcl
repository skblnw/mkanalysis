set all [atomselect top "all"] 
set some [atomselect top "resid 1 to 5"] 
measure sasa 1.4 $all -restrict $some -points sasapoints 
foreach pt $sasapoints { 
  draw point $pt 
}