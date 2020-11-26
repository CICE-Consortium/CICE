# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        gx3ncarbulk,diag1,run5day       
smoke          gx3     1x1        gx3ncarbulk,diag1,run1day       
restart        gbox128 8x1        diag1
restart        gx3     4x2        gx3ncarbulk,debug,diag1,run5day
smoke          gx3     4x1        gx3ncarbulk,diag1,run5day,thread   smoke_gx3_8x2_diag1_gx3ncarbulk_run5day

# jra55
# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day
smoke          gx3     1x1        diag1,run1day
restart        gbox128 8x1        diag1
restart        gx3     4x2        debug,diag1,run5day
smoke          gx3     4x1        diag1,run5day,thread   smoke_gx3_8x2_diag1_run5day
