i="$1"
BIN=400

g_wham_d -ip pdofilelist.dat -bins $BIN -min 6 -max 10 -temp 310 -b 0 -zprof0 9.8 -unit kCal -o pmf-$i-$BIN.xvg -hist histo-$i-$BIN.xvg
