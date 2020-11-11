# Test         Grid    PEs        Sets    BFB-compare
restart        gx3     4x2x25x29x4   dslenderX2,jra55_gx3
decomp         gx3     4x2x25x29x5   jra55_gx3
sleep 30
restart        gx3     1x1x50x58x4   jra55_gx3,droundrobin,thread restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3
restart        gx3     4x1x25x116x1  jra55_gx3,dslenderX1,thread  restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3
restart        gx3     6x2x4x29x18   jra55_gx3,dspacecurve        restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3 
restart        gx3     8x2x8x10x20   jra55_gx3,droundrobin        restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3 
restart        gx3     6x2x50x58x1   jra55_gx3,droundrobin        restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3 
restart        gx3     4x2x19x19x10  jra55_gx3,droundrobin        restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3 
restart        gx3     1x20x5x29x80  jra55_gx3,dsectrobin,short   restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3 
restart        gx3     16x2x5x10x20  jra55_gx3,drakeX2            restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3 
restart        gx3     8x2x8x10x20   jra55_gx3,droundrobin,maskhalo restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3   
restart        gx3     1x4x25x29x16  jra55_gx3,droundrobin        restart_gx3_4x2x25x29x4_dslenderX2_jra55_gx3 

