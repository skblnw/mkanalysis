set fp [open "str.dat" r]
set file_data [read $fp]
close $fp

set data [split $file_data "\n"]
foreach line $data {
  draw sphere $line radius 2
}
