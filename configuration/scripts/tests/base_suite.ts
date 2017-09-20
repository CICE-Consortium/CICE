# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day
smoke          gx3     8x2        diag24,run1year,medium
smoke          gx3     4x1        debug,diag1,run5day
smoke          gx3     8x2        debug,diag1,run5day
smoke          gx3     4x2        diag1,run5day          smoke_gx3_8x2_diag1_run5day
smoke          gx3     4x1        diag1,run5day,thread   smoke_gx3_8x2_diag1_run5day
smoke          gx3     4x1        diag1,run5day          smoke_gx3_4x1_diag1_run5day_thread
restart        gx3     8x1        none
restart        gx3     4x2        debug
restart        gx3     8x2        diag1,pondcesm
restart        gx3     8x2        diag1,pondtopo
