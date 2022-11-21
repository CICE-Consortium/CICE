# Test         Grid    PEs           Sets       BFB-compare
smoke          gx3     8x2           diag1,run5day
# decomp_suite
restart        gx3     4x2x25x29x4   dslenderX2
smoke          gx3     4x2x25x29x4   debug,run2day,dslenderX2
# reprosum_suite
smoke          gx3     4x2x25x29x4   dslenderX2,diag1,reprosum
# travis_suite
smoke          gx3     1x2           run2day
# gridsys_suite
smoke          gx3     1x1x100x116x1 reprosum,run10day
smoke          gx1     32x1x16x16x32 reprosum,run10day
smoke          gx3     1x1x100x116x1 reprosum,run10day,gridcd
smoke          gx1     32x1x16x16x32 reprosum,run10day,gridcd
smoke          gx3     1x1x100x116x1 reprosum,run10day,gridc
smoke          gx1     32x1x16x16x32 reprosum,run10day,gridc
# perf_suite
smoke          gx1     32x1x16x16x15 run2day,droundrobin
smoke          gx1     64x1x16x16x8  run2day,droundrobin,thread
