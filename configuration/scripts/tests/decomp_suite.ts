# Test         Grid    PEs        Sets    BFB-compare
restart        gx3     4x2x25x29x4   dslenderX2
restart        gx1     64x1x16x16x10 dwghtfile
restart        gbox180 16x1x6x6x60   dspacecurve,debugblocks
decomp         gx3     4x2x25x29x5   none
sleep 30
restart        gx3     1x1x50x58x4   droundrobin,thread     restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     4x1x25x116x1  dslenderX1,thread      restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     6x2x4x29x18   dspacecurve            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x2x8x10x20   droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     6x2x50x58x1   droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     5x2x33x23x4   droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     4x2x19x19x10  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     20x2x5x4x30   dsectrobin,short       restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x5x10x20  drakeX2                restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x2x8x10x20   droundrobin,maskhalo   restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x4x25x29x16  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x8x30x20x32  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     1x1x120x125x1 droundrobin,thread     restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x1x1x800  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x2x2x200  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x3x3x100  droundrobin            restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     16x2x8x8x80   dspiralcenter          restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     10x1x10x29x4  dsquarepop,thread      restart_gx3_4x2x25x29x4_dslenderX2
restart        gx3     8x1x25x29x4   drakeX2,thread         restart_gx3_4x2x25x29x4_dslenderX2

smoke        gx3     4x2x25x29x4   debug,run2day,dslenderX2
smoke        gx1     64x1x16x16x10 debug,run2day,dwghtfile
smoke        gbox180 16x1x6x6x60   debug,run2day,dspacecurve,debugblocks
sleep 30
smoke        gx3     1x1x25x58x8   debug,run2day,droundrobin,thread     smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     20x1x5x116x1  debug,run2day,dslenderX1,thread      smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     6x2x4x29x18   debug,run2day,dspacecurve            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     8x2x10x12x16  debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     6x2x50x58x1   debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     5x2x33x23x4   debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     4x2x19x19x10  debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     20x2x5x4x30   debug,run2day,dsectrobin,short       smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     16x2x5x10x20  debug,run2day,drakeX2                smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     8x2x8x10x20   debug,run2day,droundrobin,maskhalo   smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     1x6x25x29x16  debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     1x8x30x20x32  debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     1x1x120x125x1 debug,run2day,droundrobin,thread     smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     16x2x1x1x800  debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     16x2x2x2x200  debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     16x2x3x3x100  debug,run2day,droundrobin            smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     16x2x8x8x80   debug,run2day,dspiralcenter          smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     10x1x10x29x4  debug,run2day,dsquarepop,thread      smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
smoke        gx3     8x1x25x29x4   debug,run2day,drakeX2,thread         smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day

