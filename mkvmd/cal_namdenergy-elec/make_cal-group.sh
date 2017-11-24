# region index
# rr=1
# region1 
# sequence="1 3 5 7 9 11 13 15 17 19 21 23"
# sequence="5 7"
# sel1="103 107"
# sel2="150 151 154"
# pair=1
# sequence="1 3 5 7 9 11 13 15 17 19 21 23"
# sel1="150 151 154"
# sel2="103 107"
# pair=2
# sequence="1 3 5 7 9 11 13 15 17 19 21 23"
# sel1="122 125 126 129"
# sel2="122 125 126 129"
# pair=3
# sequence="2 4 6 8 10 12 14 16 18 20 22 24"
# sequence="6 8"
# sel1="247 248 244 243 240 239 236"
# sel2="236 239 240 243 244 248 247"
# pair=4

# region2
# sequence="8 4 20 8 2 6 6 18 24 16 16 12 10 14 14 22"
# sequence="6 18"
# sel1="278 280 281"
# sel2="5 6 8"
# pair=1
# sequence="8 3 20 7 2 5 6 17 24 15 16 11 10 13 14 21"
# sequence="6 17"
# sel1="277 281"
# sel2="234 237 241"
# pair=2
# sequence="8 3 20 7 2 5 6 17 24 15 16 11 10 13 14 21"
# sel1="336 335 337 330"
# sel2="240 236 237 233"
# pair=3

# region3
# sequence="16 6"
# sel1="310 311 313 335 334"
# sel2="82 81 78 85 86 89 208 209"
# pair=1
# sequence="11 6"
# sel1="236 240"
# sel2="92 96 90"
# pair=2


if true; then

    egrep -v "^#" LINES > TMP
    while read line; do
        iregion=`echo $line | awk '{print $1}'`
        ipair=`echo $line | awk '{print $2}'`
        case $iregion in
        1)
            case $ipair in
            1)
                sequence="7 3 19 7 1 5 5 17 23 15 15 11 9 13 13 21"
            ;;
            [2-3])
                sequence="8 4 20 8 2 6 6 18 24 16 16 12 10 14 14 22"
            ;;
            [4-6])
                sequence="8 3 20 7 2 5 6 17 24 15 16 11 10 13 14 21"
            ;;
            esac
            ;;
        2)
            case $ipair in
            [1-3])
                sequence="1 3 5 7 9 11 13 15 17 19 21 23";;
            4)
                sequence="2 4 6 8 10 12 14 16 18 20 22 24";;
            esac;;
        3)
            case $ipair in
            1)
                sequence="12 5 6 11 16 17 18 15";;
            [2-6])
                sequence="2 12 12 2 6 16 16 6 18 24 24 18";;
            esac
        esac

        sed -e "s/SEQ/${sequence}/g" \
            template-cal-group-tcl > cal_group.tcl
        # vmd -e <tcl> -args <psf> <dcd> <line>
        vmd -dispdev text -e cal_group.tcl -args ../combine_dcd/sel-region${iregion}.psf ../combine_dcd/mdff-region${iregion}-s50.dcd $line </dev/null
        wait
    done < TMP
    rm -f TMP

    wait
fi
