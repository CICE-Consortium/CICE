# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     1x2        diag1,run5day
smoke          gx3     2x1        debug,diag1,run5day
smoke          gx3     1x2        debug,diag1,run5day
smoke          gx3     1x1        diag1,run5day,thread   smoke_gx3_1x2_diag1_run5day
smoke          gx3     2x1        diag1,run5day,thread   smoke_gx3_1x2_diag1_run5day
restart        gx3     2x1        diag1
restart        gx3     1x2        diag1
restart        gx3     2x1        diag1,pondcesm
restart        gx3     2x1        diag1,pondtopo
