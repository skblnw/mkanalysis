# region index
ii=1
# region1 
sequence="1 3 5 7 9 11 13 15 17 19 21 23"
sel1="103 107"
sel2="150 151 154"
sequence="1 3 5 7 9 11 13 15 17 19 21 23"
sel1="150 151 154"
sel2="103 107"
sequence="1 3 5 7 9 11 13 15 17 19 21 23"
sel1="122 125 126 129"
sel2="122 125 126 129"
sequence="2 4 6 8 10 12 14 16 18 20 22 24"
sel1="247 248 244 243"
sel2="236 239 240 243 244"
sequence="2 4 6 8 10 12 14 16 18 20 22 24"
sel1="244 243 240 239 236"
sel2="243 244 248 247"

# region2
sequence="8 4 20 8 2 6 6 18 24 16 16 12 10 14 14 22"
sel1="278 280 281"
sel2="5 6 8"
sequence="8 3 20 7 2 5 6 17 24 15 16 11 10 13 14 21"
sel1="277 281"
sel2="234 237 241"
sequence="8 3 20 7 2 5 6 17 24 15 16 11 10 13 14 21"
sel1="336 335 337 330"
sel2="240 236 237 233"


if false; then
    sed -e 's/SEQ/'$sequence'/g' \
        -e 's/TEXT1/'$sel1'/g' \
        -e 's/TEXT2/'$sel2'/g' \
        template-cal-energy-tcl > cal_energy.tcl

    # vmd -e cal_energy.tcl -args <psf> <dcd> <outname>
    vmd -e cal_energy.tcl -args ../combine_dcd/sel-region$ii.psf ../combine_dcd/mdff-region$ii-s50.dcd region$ii
    wait
fi

# region 3
if false; then
sed -e 's/SEQ/2 12 16 6 6 16 24 18/g' \
    -e 's/TEXT1/236 239 240/g' \
    -e 's/TEXT1/90 92 96/g' \
    template-cal-energy-tcl > cal_energy.tcl
fi