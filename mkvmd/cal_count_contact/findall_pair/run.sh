seg1="4 10 1 3 8 12 5 7"
seg2="2 4 3 9 6 8 7 11"
sed -e 's/SEG1/'"$seg1"'/g' -e 's/SEG2/'"$seg2"'/g' template-measure-findall-tcl > vm_measure-findall-1.tcl

seg1="1 3 5 7 9 11"
seg2="2 4 6 8 10 12"
sed -e 's/SEG1/'"$seg1"'/g' -e 's/SEG2/'"$seg2"'/g' template-measure-findall-tcl > vm_measure-findall-2.tcl

seg1="1 6 3 8 9"
seg2="6 3 8 9 12"
sed -e 's/SEG1/'"$seg1"'/g' -e 's/SEG2/'"$seg2"'/g' template-measure-findall-tcl > vm_measure-findall-3.tcl

for ii in {1..3}
do
    sed -e 's/NAME/'$ii'/g' template-run-sh > run-$ii.sh
    nohup sh run-$ii.sh >& output-md/LOG-run-$ii.log &
done
