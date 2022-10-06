# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day,e3sm
smoke          gx3     1x1        diag1,run1day,e3smbgc
restart        gbox128 8x1        diag1,e3sm
restart        gx3     4x2        debug,diag1,e3smbgc
smoke          gx3     4x1        diag1,run5day,thread,e3sm   smoke_gx3_8x2_diag1_run5day
