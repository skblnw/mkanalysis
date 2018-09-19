#!/bin/bash

# VMD binary path
# version > 1.9.2
# TachyonLOptiXInternal compiled
VMD=vmd

WRAPPER_NAME=rwrapper
OUTNAME=tachyon

# Create a render wrapper
cat render.vmd > $WRAPPER_NAME

# Add render scripts
cat >> $WRAPPER_NAME << EOF

light 0 on
light 1 off
light 2 off
light 3 off
display resize 953 964
display shadows on
display ambientocclusion on
display aoambient 0.60
display aodirect 0.60

render aasamples TachyonLOptiXInternal 12
render TachyonLOptiXInternal $OUTNAME.tga
# Convert tga to jpg using ImageMagick convert
# You may put it in the mother BASH instead
#exec convert $OUTNAME.tga $DIR_OUTPUT/$OUTNAME.jpg

EOF

# Main rendering process
$VMD -dispdev text -eofexit < $WRAPPER_NAME
wait
