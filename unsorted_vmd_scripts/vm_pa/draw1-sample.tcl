    # draw some nice vectors
    #graphics $mol delete all
    graphics $mol color yellow
    graphics $mol material EdgyShiny
    set COM [Orient::sel_com $sel $weights]
    vmd_draw_vector $mol $COM [vecscale $scale $a1]
    vmd_draw_vector $mol $COM [vecscale [vecscale $scale $a2] -1]
    vmd_draw_vector $mol $COM [vecscale $scale $a3]

    graphics $mol color black
    graphics $mol text [vecadd $COM [vecscale $scale2 $a1]] "1" size 3 thickness 6
    graphics $mol text [vecadd [vecadd $COM [vecscale [vecscale $scale2 $a2] -1]] {5 0 -8}] "2" size 3 thickness 6
    graphics $mol text [vecadd [vecadd $COM [vecscale $scale2 $a3]] {-20 -20 10}] "3" size 3 thickness 6