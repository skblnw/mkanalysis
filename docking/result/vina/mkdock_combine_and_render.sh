#!/bin/bash

DIR_DOCKRES=
RECEPTOR=
DIR_RENDER=
DIR_OUTPUT=

for filename in $(cat namelist)
do
    ZINC=$(basename $filename)
    echo $ZINC
    mkdir -p $ZINC
    cd $ZINC
    sed -n '/^MODEL 1$/,/^ENDMDL$/p' $DIR_DOCKRES/$ZINC/out.pdbqt > tmp
    grep '^ATOM' tmp > tmp2
    mv tmp2 tmp

    head -n -1 $RECEPTOR > combine.pdb
    cat tmp >> combine.pdb
    echo -e "TER\nENDMDL" >> combine.pdb

    cp $DIR_RENDER/render_apo.vmd render.vmd
    cp $DIR_RENDER/make_render_vmd.sh make.sh
    bash make.sh
    convert tachyon.tga $DIR_OUTPUT/${ZINC}_apo.jpg
    mv combine.pdb combine_apo.pdb
    rm -f tmp render.vmd make.sh

    cd ..
done