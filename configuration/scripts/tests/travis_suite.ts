# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     1x2        zzgx3ncarbulk,run2day
smoke          gx3     1x1        zzgx3ncarbulk,debug,run1day
smoke          gx3     2x2        zzgx3ncarbulk,debug,run1day
smoke          gx3     2x1        zzgx3ncarbulk,run2day,thread  smoke_gx3_1x2_run2day_zzgx3ncarbulk
restart        gx3     2x1        zzgx3ncarbulk
restart        gx3     1x2        zzgx3ncarbulk

# jra55
# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     1x2        run2day
smoke          gx3     1x1        debug,run1day
smoke          gx3     2x2        debug,run1day
smoke          gx3     2x1        run2day,thread  smoke_gx3_1x2_run2day
restart        gx3     2x1        
restart        gx3     1x2        
