# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     2x2        debug,diag1,run1day
smoke          gx3     1x2        debug,diag1,run1day smoke_gx3_2x2_debug_diag1_run1day
smoke          gx3     1x1        debug,diag1,run1day,thread smoke_gx3_2x2_debug_diag1_run1day
smoke          gx3     1x1        diag1,run5day
smoke          gx3     2x1        diag1,run5day smoke_gx3_1x1_diag1_run5day
restart        gx3     1x1        diag1
restart        gx3     2x2        diag1
