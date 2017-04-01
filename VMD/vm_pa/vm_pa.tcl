# 1. source vm_pa.tcl to do all the things
# 2. configure draw1-sample.tcl for arrows of BAR and 2 PHs
# 3. orient.tcl has been modified and will be sourced on the first line of vm_pa.tcl
source orient.tcl

set selall [atomselect top all]
set selbar [atomselect top "resid 1 to 250 and name CA"]
draw delete all
set I1 [draw principalaxes $selbar 1]
set selph1 [atomselect top "segname P5 and resid 345 to 361 and name CA"]
set I2 [draw principalaxes $selph1 2]
set selph2 [atomselect top "segname P6 and resid 345 to 361 and name CA"]
set I3 [draw principalaxes $selph2 3]