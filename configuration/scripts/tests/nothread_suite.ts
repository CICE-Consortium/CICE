# Test         Grid    PEs        Sets    BFB-compare

restart        gx3     8x1x25x29x2   dslenderX2
logbfb         gx3     8x1x25x29x2   dslenderX2,diag1,reprosum

smoke          gx3     16x1       diag1,run5day
smoke          gx3     1x1        debug,diag1,run2day
smoke          gx3     4x1        debug,diag1,run5day,thread
restart        gx3     16x1       thread
smoke          gx3     16x1       diag24,run1year,medium

#tcraig, hangs nodes intermittently on izumi
#restart        gx1     160x1      droundrobin,medium
#restart        tx1     160x1      dsectrobin,medium

restart        gx3     16x1       none
restart        gx3     16x1       gx3ncarbulk,iobinary

restart        gx3     12x1       alt01
restart        gx3     16x1       alt02
restart        gx3     8x1        alt03
restart        gx3     16x1       alt04
restart        gx3     16x1       alt05
restart        gx3     20x1       alt06
restart        gx3     18x1       alt01,debug,short
restart        gx3     20x1       alt02,debug,short
restart        gx3     24x1       alt03,debug,short
smoke          gx3     24x1       alt04,debug,short
smoke          gx3     32x1       alt05,debug,short
smoke          gx3     16x1       alt06,debug,short
restart        gx3     16x1       isotope
smoke          gx3     6x1        isotope,debug
smoke          gx3     8x1        fsd1,diag24,run5day,debug
smoke          gx3     16x1       fsd12,diag24,run5day,short
restart        gx3     12x1       fsd12,debug,short
smoke          gx3     20x1       fsd12ww3,diag24,run1day,medium

restart        gbox128 8x1        short
restart        gbox128 16x1       boxnodyn,short
restart        gbox128 24x1       boxnodyn,short,debug
restart        gbox128 12x1       boxadv,short
smoke          gbox128 20x1       boxadv,short,debug
restart        gbox128 32x1       boxrestore,short
smoke          gbox128 24x1       boxrestore,short,debug
restart        gbox80  1x1        box2001
smoke          gbox80  1x1        boxslotcyl

smoke          gx3     16x1        medium,run90day,yi2008
restart        gx3     12x1        short
#tcraig, hangs nodes intermittently on izumi
#smoke          gx1     24x1       medium,run90day,yi2008
#restart        gx1     24x1       short

smoke          gx3     16x1       bgcz
smoke          gx3     16x1       bgcz,debug
smoke          gx3     24x1       bgcskl,debug
#tcraig, hangs nodes intermittently on izumi
#restart        gx1     128x1      bgcsklclim,medium
#restart        gx1     256x1      bgczclim,medium

decomp         gx3     8x1x5x29x20   none
restart        gx3     1x1x50x58x4   droundrobin        restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     4x1x25x116x1  dslenderX1         restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     12x1x4x29x9   dspacecurve        restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     16x1x8x10x10  droundrobin        restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     6x1x50x58x1   droundrobin        restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     8x1x19x19x5   droundrobin        restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     20x1x5x29x20  dsectrobin,short   restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     32x1x5x10x10  drakeX2            restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     16x1x8x10x10  droundrobin,maskhalo   restart_gx3_8x1x25x29x2_dslenderX2
restart        gx3     4x1x25x29x4   droundrobin        restart_gx3_8x1x25x29x2_dslenderX2

logbfb         gx3     1x1x50x58x4   droundrobin,diag1,maskhalo,reprosum logbfb_gx3_8x1x25x29x2_diag1_dslenderX2_reprosum
logbfb         gx3     4x1x25x116x1  dslenderX1,diag1,maskhalo,reprosum  logbfb_gx3_8x1x25x29x2_diag1_dslenderX2_reprosum
logbfb         gx3     20x1x5x29x20  dsectrobin,diag1,short,reprosum     logbfb_gx3_8x1x25x29x2_diag1_dslenderX2_reprosum
logbfb         gx3     16x1x8x10x10  droundrobin,diag1,reprosum          logbfb_gx3_8x1x25x29x2_diag1_dslenderX2_reprosum
logbfb         gx3     6x1x50x58x1   droundrobin,diag1,reprosum          logbfb_gx3_8x1x25x29x2_diag1_dslenderX2_reprosum
logbfb         gx3     12x1x4x29x9   dspacecurve,diag1,maskhalo,reprosum logbfb_gx3_8x1x25x29x2_diag1_dslenderX2_reprosum
