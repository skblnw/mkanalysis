#!/bin/bash

DIR_DOCKRES=output_aa
RECEPTOR=receptor.pdbqt

ls $DIR_DOCKRES > namelist
for filename in $(cat namelist)
do
    ZINC=$(basename $filename)
    sed -n '/^MODEL 1$/,/^ENDMDL$/p' $DIR_DOCKRES/$ZINC/out.pdbqt > tmp
    grep '^HETATM' tmp > tmp2
    mv tmp2 tmp

    head -n -1 $RECEPTOR > combine.pdb
    cat tmp >> combine.pdb
    echo -e "TER\nENDMDL" >> combine.pdb
    mv combine.pdb $DIR_DOCKRES/$ZINC
done
rm -f tmp
