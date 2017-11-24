DCDIN=cls1-gs3-water-protein-35ns-s50.dcd
LINES1=/home/kevin/Dropbox/QWD/BAR-PH/mdff/run-cls1-gs3-water/run/make_analysis/cal_contact_pairs/LINES-I1
LINES2=/home/kevin/Dropbox/QWD/BAR-PH/mdff/run-cls1-gs3-water/run/make_analysis/cal_contact_pairs/LINES-I2

nn=0
IFS=$'\n'
grep --no-filename -v "^#" $LINES1 | awk '{print "segname P1 P2 P5 P6 and resid "$1" "$2" "$3" or segname P3 P4 P7 P8 and resid "$4" "$5" "$6}' | while read -r SELTEXT
do
((nn++))
OUTNAME=I1-sel$nn
IDXOUT=index-$OUTNAME
PDBOUT=initial-$OUTNAME
DCDOUT=output/`sed 's/protein/'$OUTNAME'/g' <<< $DCDIN`

sed -e 's/SELTEXT/'$SELTEXT'/g' -e 's/IDXOUT/'$IDXOUT'/g' -e 's/PDBOUT/'$PDBOUT'/g' template-get_index_and_pdb > vm_get.tcl

vmd -dispdev text -e vm_get.tcl < /dev/null
wait

catdcd -o $DCDOUT -i output/$IDXOUT.txt $DCDIN
wait
done

grep --no-filename -v "^#" $LINES2 | awk '{print "segname P1 P2 P7 P8 and resid "$1" "$2" "$3" or segname P5 P6 P3 P4 and resid "$4" "$5" "$6}' | while read -r SELTEXT
do
((nn++))
OUTNAME=I2-sel$nn
IDXOUT=index-$OUTNAME
PDBOUT=initial-$OUTNAME
DCDOUT=output/`sed 's/protein/'$OUTNAME'/g' <<< $DCDIN`

sed -e 's/SELTEXT/'$SELTEXT'/g' -e 's/IDXOUT/'$IDXOUT'/g' -e 's/PDBOUT/'$PDBOUT'/g' template-get_index_and_pdb > vm_get.tcl

vmd -dispdev text -e vm_get.tcl < /dev/null
wait

catdcd -o $DCDOUT -i output/$IDXOUT.txt $DCDIN
wait
done
