# Test         Grid    PEs        Sets    BFB-compare
restart        gx3     4x2x25x29x4   dslenderX2
restart        gx1     64x1x16x16x10 dwghtfile
restart        gx1     32x2x10x12x32 dsectcart,short
restart        gbox180 16x1x6x6x60   dspacecurve,debugblocks
restart        gbox80  4x2x23x21x6   boxgauss,bczerogradient
restart        gbox80  2x2x29x29     boxgauss,bclinearextrap
restart        gbox80  3x2x22x21     boxgauss,bccyclicextrap
restart        gbox80  6x2x13x12     boxgauss,bcclosed
restart        gbox80  8x2x6x7       boxgauss,bccyclic
decomp         gx3     4x2x25x29x5   none
decomp         gx3     4x2x25x29     none
decomp         gx3     4x2x25x29x5   dynpicard,reprosum
decomp         gx3     4x2x25x29x5   dyneap
restart        gx3     1x1x50x58x4   droundrobin,thread     restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     4x1x25x116x1  dslenderX1,thread      restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     6x2x4x29x18   dspacecurve            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x2x8x10x20   droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     6x2x50x58x1   droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     5x2x33x23x4   droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     4x2x19x19x10  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     20x2x5x4x30   dsectrobin,short       restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x5x10     drakeX2                restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x2x8x10x20   droundrobin,maskhalo   restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x4x25x29x16  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x8x30x20x32  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x1x120x125x1 droundrobin,thread     restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x1x1x800  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x2x2x200  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x3x3x100  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x8x8x80   dspiralcenter          restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     10x1x10x29x4  dsquarepop,thread      restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x1x25x29     drakeX2,thread         restart_gx3_4x2x25x29x4_dslenderX2
restart        gbox80  1x2x23x21x20  boxgauss,bczerogradient         restart_gbox80_4x2x23x21x6_bczerogradient_boxgauss
restart        gbox80  1x1x15x17     boxgauss,bclinearextrap,thread  restart_gbox80_2x2x29x29_bclinearextrap_boxgauss
restart        gbox80  1x2x29x28     boxgauss,bccyclicextrap         restart_gbox80_3x2x22x21_bccyclicextrap_boxgauss
restart        gbox80  1x4x12x13     boxgauss,bcclosed          restart_gbox80_6x2x13x12_bcclosed_boxgauss
restart        gbox80  1x8x9x8       boxgauss,bccyclic          restart_gbox80_8x2x6x7_bccyclic_boxgauss

smoke        gx3     4x2x25x29     debug,run2day,dslenderX2
smoke        gx1     64x1x16x16    debug,run2day,dwghtfile
smoke        gx1     32x2x10x12    debug,run2day,dsectcart
smoke        gbox180 16x1x6x6      debug,run2day,dspacecurve,debugblocks
smoke        gbox80  4x2x23x21x6   debug,run2day,boxgauss,bczerogradient
smoke        gbox80  3x2x22x21     debug,run2day,boxgauss,bccyclicextrap
smoke        gbox80  2x2x29x29     debug,run2day,boxgauss,bclinearextrap
smoke        gx3     1x1x25x58     debug,run2day,droundrobin,thread     smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     20x1x5x116    debug,run2day,dslenderX1,thread      smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     6x2x4x29      debug,run2day,dspacecurve            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     8x2x10x12x18  debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     6x2x50x58     debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     5x2x33x23     debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     4x2x19x19x10  debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     20x2x5x4      debug,run2day,dsectrobin,short       smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     16x2x5x10     debug,run2day,drakeX2                smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     8x2x8x10x20   debug,run2day,droundrobin,maskhalo   smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     1x6x25x29x16  debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     1x8x30x20x32  debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     1x1x120x125x1 debug,run2day,droundrobin,thread     smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     16x2x1x1x800  debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     16x2x2x2x200  debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     16x2x3x3x100  debug,run2day,droundrobin            smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     16x2x8x8      debug,run2day,dspiralcenter          smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     10x1x10x29    debug,run2day,dsquarepop,thread      smoke_gx3_4x2x25x29_debug_dslenderX2_run2day
smoke        gx3     8x1x25x29     debug,run2day,drakeX2,thread         smoke_gx3_4x2x25x29_debug_dslenderX2_run2day

