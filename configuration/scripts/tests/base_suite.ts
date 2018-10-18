# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day
smoke          gx3     1x1        debug,diag1,run2day
smoke          gx3     1x4        debug,diag1,run2day
smoke          gx3     4x1        debug,diag1,run5day
restart        gx3     8x2        debug
smoke          gx3     8x2        diag24,run1year,medium
decomp         gx3     4x2x25x29x5
smoke          gx3     4x2        diag1,run5day          smoke_gx3_8x2_diag1_run5day
smoke          gx3     4x1        diag1,run5day,thread   smoke_gx3_8x2_diag1_run5day
restart        gx1     40x4       droundrobin,short
restart        tx1     40x4       dsectrobin,short
restart        gx3     4x4        none
restart        gx3     4x4        iobinary
restart        gx3     6x2        alt01
restart        gx3     8x2        alt02
restart        gx3     4x2        alt03
restart        gx3     4x4        alt04
restart        gx3     4x4        alt05
restart        gbox128  4x2       none
restart        gbox128  4x2       boxdyn
restart        gbox128  2x2       boxadv
restart        gbox128  4x4       boxrestore
smoke          gx3     8x2        bgcISPOL,short
smoke          gx3     8x2        bgcNICE,debug
smoke          gx3     8x1        bgcNICE
smoke          gx3     4x2        bgcNICE,thread   smoke_gx3_8x1_bgcNICE
smoke          gx3     1x1        bgcNICE          smoke_gx3_4x2_bgcNICE_thread
restart        gx1     4x2        sklgx1,medium
restart        gx1     8x1        bgcgx1,medium
