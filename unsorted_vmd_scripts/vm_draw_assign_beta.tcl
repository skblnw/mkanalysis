# r: starting residue #
set r 1
foreach i {0.8	2.9	4.6	5.4	5.7	4.6	2.6	2.9	2.1	0.0	0.0	0.0	0.0	0.0	0.0	2.5	1.6	0.0	0.3	2.2	0.1	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	} {
    set sel [atomselect top "segname P1 and resid ${r}"]
	$sel set beta ${i}
	incr r
}