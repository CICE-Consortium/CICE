# Test         Grid    PEs           Sets       BFB-compare
smoke          gx3     8x2           diag1,run5day
restart        gx3     4x2x25x29x4   dslenderX2
smoke          gx3     4x2x25x29x4   debug,run2day,dslenderX2
smoke          gx3     4x2x25x29x4   dslenderX2,diag1,reprosum,cmplog
smoke          gx3     1x2           run2day
smoke          gx3     1x1x100x116x1   reprosum,run10day
smoke          gx1     32x1x16x16x32   reprosum,run10day
smoke          gx3     1x1x100x116x1   reprosum,run10day,gridcd
smoke          gx1     32x1x16x16x32   reprosum,run10day,gridcd
smoke          gx3     1x1x100x116x1   reprosum,run10day,gridc
smoke          gx1     32x1x16x16x32   reprosum,run10day,gridc
