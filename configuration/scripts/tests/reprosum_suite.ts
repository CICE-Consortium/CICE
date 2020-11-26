# Test         Grid    PEs        Sets    BFB-compare
logbfb        gx3     4x2x25x29x4   zzgx3ncarbulk,dslenderX2,diag1,reprosum
#logbfb        gx3     4x2x25x29x4   zzgx3ncarbulk,dslenderX2,diag1
sleep 60
logbfb        gx3     1x1x50x58x4   zzgx3ncarbulk,droundrobin,diag1,thread,maskhalo,reprosum logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum_zzgx3ncarbulk
logbfb        gx3     4x1x25x116x1  zzgx3ncarbulk,dslenderX1,diag1,thread,maskhalo,reprosum  logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum_zzgx3ncarbulk
logbfb        gx3     1x20x5x29x80  zzgx3ncarbulk,dsectrobin,diag1,short,reprosum   logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum_zzgx3ncarbulk
logbfb        gx3     8x2x8x10x20   zzgx3ncarbulk,droundrobin,diag1,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum_zzgx3ncarbulk
logbfb        gx3     6x2x50x58x1   zzgx3ncarbulk,droundrobin,diag1,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum_zzgx3ncarbulk
logbfb        gx3     6x2x4x29x18   zzgx3ncarbulk,dspacecurve,diag1,maskhalo,reprosum        logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_reprosum_zzgx3ncarbulk
#logbfb        gx3     8x2x8x10x20   zzgx3ncarbulk,droundrobin,diag1                 logbfb_gx3_4x2x25x29x4_diag1_dslenderX2_zzgx3ncarbulk

# jra55
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
