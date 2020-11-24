# Test         Grid    PEs        Sets    BFB-compare
restart        gx3     16x1       diag1
smoke          gx3     1x1        debug,diag1,run2day
smoke          gx3     4x1        debug,diag1,run2day,thread

# jra55
# Test         Grid    PEs        Sets    BFB-compare
restart        gx3     16x1       jra55_gx3,diag1
smoke          gx3     1x1        jra55_gx3,debug,diag1,run2day
smoke          gx3     4x1        jra55_gx3,debug,diag1,run2day,thread

