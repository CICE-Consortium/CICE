# Test         Grid    PEs           Sets       BFB-compare
smoke          gx3     8x2           diag1,run5day
restart        gx3     4x2x25x29x4   dslenderX2
logbfb         gx3     4x2x25x29x4   dslenderX2,diag1,reprosum
smoke          gx3     1x2           run2day
