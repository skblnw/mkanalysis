for ii in 616 630 640 650 660 670
#for ii in 680 695 710 725 740 755 765 770 785
#for ii in 800 815 830 845 860 875 890 905 920 935 950 965 980
do
./catdcd -o output/us-z$ii-s100.dcd -i index.txt -stride 100 ../../output/us-z$ii-{13..21}.dcd
done
