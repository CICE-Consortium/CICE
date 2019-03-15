# Test         Grid    PEs        Sets    BFB-compare
restart        gx3     4x2x25x29x4   dslenderX2,reprosum
decomp         gx3     4x2x25x29x5   reprosum
restart        gx3     4x1x25x116x1  dslenderX1,thread,reprosum  restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     6x2x4x29x18   dspacecurve,reprosum        restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     8x2x8x10x20   droundrobin,reprosum        restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     6x2x50x58x1   droundrobin,reprosum        restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     4x2x19x19x10  droundrobin,reprosum        restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     1x20x5x29x80  dsectrobin,short,reprosum   restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     16x2x5x10x20  drakeX2,reprosum            restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     8x2x8x10x20   droundrobin,maskhalo,reprosum   restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     1x4x25x29x16  droundrobin,reprosum        restart_gx3_4x2x25x29x4_dslenderX2_reprosum
restart        gx3     1x1x50x58x4   droundrobin,thread,reprosum restart_gx3_4x2x25x29x4_dslenderX2_reprosum
