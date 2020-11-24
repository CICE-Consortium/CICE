# Test         Grid    PEs        Sets    BFB-compare
logbfb        gx3     4x2x25x29x4   dslenderX2,diag1,reprosum
#logbfb        gx3     4x2x25x29x4   dslenderX2,diag1
sleep 60
logbfb        gx3     1x1x50x58x4   droundrobin,diag1,thread,maskhalo,reprosum logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum
logbfb        gx3     4x1x25x116x1  dslenderX1,diag1,thread,maskhalo,reprosum  logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum
logbfb        gx3     1x20x5x29x80  dsectrobin,diag1,short,reprosum   logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum
logbfb        gx3     8x2x8x10x20   droundrobin,diag1,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum
logbfb        gx3     6x2x50x58x1   droundrobin,diag1,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum
logbfb        gx3     6x2x4x29x18   dspacecurve,diag1,maskhalo,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum
#logbfb        gx3     8x2x8x10x20   droundrobin,diag1                 logbfb_gx3_4x2x25x29x4_diag1_dslenderX2

# jra55
# Test         Grid    PEs        Sets    BFB-compare
logbfb        gx3     4x2x25x29x4   jra55_gx3,dslenderX2,diag1,reprosum
#logbfb        gx3     4x2x25x29x4   jra55_gx3,dslenderX2,diag1
sleep 60
logbfb        gx3     1x1x50x58x4   jra55_gx3,droundrobin,diag1,thread,maskhalo,reprosum logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_jra55_gx3_reprosum
logbfb        gx3     4x1x25x116x1  jra55_gx3,dslenderX1,diag1,thread,maskhalo,reprosum  logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_jra55_gx3_reprosum
logbfb        gx3     1x20x5x29x80  jra55_gx3,dsectrobin,diag1,short,reprosum   logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_jra55_gx3_reprosum
logbfb        gx3     8x2x8x10x20   jra55_gx3,droundrobin,diag1,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_jra55_gx3_reprosum
logbfb        gx3     6x2x50x58x1   jra55_gx3,droundrobin,diag1,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_jra55_gx3_reprosum
logbfb        gx3     6x2x4x29x18   jra55_gx3,dspacecurve,diag1,maskhalo,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_jra55_gx3_reprosum
#logbfb        gx3     8x2x8x10x20   jra55_gx3,droundrobin,diag1                 logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_jra55_gx3
