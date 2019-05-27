#!/bin/bash

# VMD binary path
# version > 1.9.2
# TachyonLOptiXInternal compiled
VMD=vmd

# VS: visualisation state
# Where you put your DCD and VS files
DCD_DIR=/home/chan773/playground/cpf1/run/output
# Where you wanna output JPG
RENDER_DIR=/home/chan773/playground/render/cpf1-amber-test
# Temporary folder for TGA files, usually leave as defualt
TMP_DIR=/home/PHARMACY/chan.773/playground/render/tmp

cd $RENDER_DIR
for name in front
do
    # Name of VS
    VS_NAME=$DCD_DIR/vs_render_$name

    OUTPUT_DIR1=$TMP_DIR/output-$name
    OUTPUT_DIR2=$RENDER_DIR/output-$name
    WRAPPER_NAME=$RENDER_DIR/rwrapper-$name

    # Create output directories
    mkdir -p $OUTPUT_DIR1 $OUTPUT_DIR2

    # Create a render wrapper
    cat $VS_NAME.vmd > $WRAPPER_NAME

    # Add frame counter
    cat >> $WRAPPER_NAME << EOF
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
    set nsperframe 1
    set time [format "%3.2f ns" [expr \$vmd_frame([molinfo top]) * \$nsperframe]] 
    draw text {-100 -60 0} "\$time" size 1 thickness 2
}
#enabletrace
EOF
	
    # Add render scripts
    cat >> $WRAPPER_NAME << EOF
axes location off
light 0 on
light 1 off
light 2 off
light 3 off
display shadows on
display resize 1920 1080
display ambientocclusion on
display aoambient 0.50
display aodirect 0.50
display update

set outdir1 $OUTPUT_DIR1
set outdir2 $OUTPUT_DIR2
set nf [molinfo top get numframes]
for {set ii 18} {\$ii < \$nf} {incr ii 1} {
        display update
        animate goto \$ii
        set filename snap_[format "%04d" [expr \$ii]]
        render aasamples TachyonLOptiXInternal 12
        render TachyonLOptiXInternal \$outdir1/\$filename.tga
        #render aasamples TachyonLOptiXInternal 12
        #render TachyonInternal \$outdir1/\$filename.tga
        puts "From \$outdir1/\$filename.tga"
        puts "To \$outdir2/\$filename.jpg"
        # Convert tga to jpg using ImageMagick convert
        exec convert \$outdir1/\$filename.tga \$outdir2/\$filename.jpg
}
EOF

    #dos2unix $WRAPPER_NAME

    # Main rendering process
    cd $DCD_DIR
    $VMD -dispdev text -eofexit < $WRAPPER_NAME
    wait

done
