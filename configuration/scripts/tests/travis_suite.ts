# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     1x2        gx3ncarbulk,run2day
smoke          gx3     1x1        gx3ncarbulk,debug,run1day
smoke          gx3     2x2        gx3ncarbulk,debug,run1day
smoke          gx3     2x1        gx3ncarbulk,run2day,thread  smoke_gx3_1x2_gx3ncarbulk_run2day
restart        gx3     2x1        gx3ncarbulk
restart        gx3     1x2        gx3ncarbulk

# jra55
# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     1x2        run2day
smoke          gx3     1x1        debug,run1day
smoke          gx3     2x2        debug,run1day
smoke          gx3     2x1        run2day,thread  smoke_gx3_1x2_run2day
restart        gx3     2x1        
restart        gx3     1x2        
