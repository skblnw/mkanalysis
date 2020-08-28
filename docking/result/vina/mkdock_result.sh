#!/bin/bash

DIR_DOCKRES=output_aa

for name in `ls $DIR_DOCKRES`
do
    echo "> grepping results for lignad: $name"
    ener=`grep "^   1 " $DIR_DOCKRES/$name/log.txt | awk {'print $2'}`
    echo $name $ener >> tmp
done

sort -h -k2 tmp > ${DIR_DOCKRES}_sort
mv tmp ${DIR_DOCKRES}_unsort
rm -f tmp