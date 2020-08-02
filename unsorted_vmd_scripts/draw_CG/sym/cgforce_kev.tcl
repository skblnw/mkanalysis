#########################################
## A short script to calculate forces from namd trajectory dcd file
## Jun Oct 2013
## Modified by Kevin Apr 2014 to output NAMD force
## Usage: (1) copy all the parameters from namd production run to the new config file
##        (2) paste the following to the new config file
##        (3) modify input/output file names, number of atoms, number of frames
##        (4) run namd with the new config file
## Output: atom number, X(3), F(3). 
## Units:  X Ang, F Kcal/Mol Ang^{-1}
## Other Notes: no virial, no velocity
## Other Notes: associate with NPT.namd
## Other Notes: unsym
#########################################

set stt   1
set natms 308201

for {set i $stt }  {$i<= $natms} {incr i} {
addatom $i
}
 
set bdindexs(1) 1
set bdindexs(2) 401
set bdindexs(3) 491
set bdindexs(4) 880
set bdindexs(5) 1615
set bdindexs(6) 2109
set bdindexs(7) 2373
set bdindexs(8) 2559
set bdindexs(9) 2724
set bdindexs(10) 3046
set bdindexs(11) 3467
set bdindexs(12) 3907
set bdindexs(13) 4274
set bdindexs(14) 4509
set bdindexs(15) 4767
set bdindexs(16) 4922
set bdindexs(17) 5091
set bdindexs(18) 5264
set bdindexs(19) 5461
set bdindexs(20) 5582
set bdindexs(21) 5904
set bdindexs(22) 6314
set bdindexs(23) 6396
set bdindexs(24) 6793
set bdindexs(25) 7519
set bdindexs(26) 8013
set bdindexs(27) 8294
set bdindexs(28) 8470
set bdindexs(29) 8628
set bdindexs(30) 8961
set bdindexs(31) 9392
set bdindexs(32) 9827
set bdindexs(33) 10144
set bdindexs(34) 10372
set bdindexs(35) 10616
set bdindexs(36) 10807
set bdindexs(37) 10954
set bdindexs(38) 11142
set bdindexs(39) 11313
set bdindexs(40) 11445

set bdindexe(1) 400
set bdindexe(2) 490
set bdindexe(3) 879
set bdindexe(4) 1614
set bdindexe(5) 2108
set bdindexe(6) 2372
set bdindexe(7) 2558
set bdindexe(8) 2723
set bdindexe(9) 3045
set bdindexe(10) 3466
set bdindexe(11) 3906
set bdindexe(12) 4273
set bdindexe(13) 4508
set bdindexe(14) 4766
set bdindexe(15) 4921
set bdindexe(16) 5090
set bdindexe(17) 5263
set bdindexe(18) 5460
set bdindexe(19) 5581
set bdindexe(20) 5903
set bdindexe(21) 6313
set bdindexe(22) 6395
set bdindexe(23) 6792
set bdindexe(24) 7518
set bdindexe(25) 8012
set bdindexe(26) 8293
set bdindexe(27) 8469
set bdindexe(28) 8627
set bdindexe(29) 8960
set bdindexe(30) 9391
set bdindexe(31) 9826
set bdindexe(32) 10143
set bdindexe(33) 10371
set bdindexe(34) 10615
set bdindexe(35) 10806
set bdindexe(36) 10953
set bdindexe(37) 11141
set bdindexe(38) 11312
set bdindexe(39) 11444
set bdindexe(40) 11725
#print $bdindexs(3) $bdindexe(4)

proc calcforces {} {
 
	global stt natms bdindexs bdindexe

	loadtotalforces f
	loadcoords coordxyz
	loadmasses matom

	if {[array exists f]} {
	  
		for {set k 1} {$k<=1} {incr k} {
			for {set j 1} {$j<=40} {incr j} {

				set summ 0
				set comx 0
				set comy 0
				set comz 0
				set tforx 0
				set tfory 0
				set tforz 0

				#set stt2 [expr {$stt1 + ($j-1)*10} ]
				#set ent2 [expr {$ent1 + ($j-1)*10} ]

				set stt2 [expr {($k-1)*6529+$bdindexs($j)} ]
				set ent2 [expr {($k-1)*6529+$bdindexe($j)} ]

				#print $k $j $stt2 $ent2

				for {set i $stt2 } {$i <= $ent2 } {incr i} {

					set m0 $matom($i)
					set f1 $f($i)
					set f1x [lindex $f1 0]
					set f1y [lindex $f1 1]
					set f1z [lindex $f1 2]
					set xyz $coordxyz($i)
					set xxx [lindex $xyz 0]
					set yyy [lindex $xyz 1]
					set zzz [lindex $xyz 2]
					set summ [expr {$summ+$m0}]
					set comx [expr {$comx+$m0*$xxx}]
					set comy [expr {$comy+$m0*$yyy}]
					set comz [expr {$comz+$m0*$zzz}]
					set tforx [expr {$tforx+$f1x}]
					set tfory [expr {$tfory+$f1y}]
					set tforz [expr {$tforz+$f1z}]
					#print $i $m0 $xxx $yyy $zzz $f1x $f1y $f1z
				 }

				set comx [expr {$comx/$summ}]
				set comy [expr {$comy/$summ}]
				set comz [expr {$comz/$summ}]

				#print "total mass, the com, and total forces are:"
				#print  $summ $comx $comy $comz $tforx $tfory $tforz
				print $comx $comy $comz $tforx $tfory $tforz
			}
		}

	}
}

