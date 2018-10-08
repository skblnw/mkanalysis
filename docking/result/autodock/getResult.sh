echo "" > result
for filename in `cat namelist`
do
    energy=`grep "^USER    Estimated Free Energy of Binding" $filename | head -1 | awk '{print $8}'`
    echo $filename $energy >> result
done
