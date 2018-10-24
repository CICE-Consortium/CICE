# Test         Grid    PEs        Sets    BFB-compare
restart        gx3     4x2x25x29x4   dslenderX2
decomp         gx3     4x2x25x29x5
restart        gx3     4x1x25x116x1  dslenderX1,thread  restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     6x2x4x29x18   dspacecurve        restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x2x8x10x20   droundrobin        restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     6x2x50x58x1   droundrobin        restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     4x2x19x19x10  droundrobin        restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x20x5x29x80  dsectrobin,short   restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x5x10x20  drakeX2            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x2x8x10x20   droundrobin,maskhalo   restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x4x25x29x16  droundrobin        restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x1x50x58x4   droundrobin,thread restart_gx3_4x2x25x29x4_dslenderX2

