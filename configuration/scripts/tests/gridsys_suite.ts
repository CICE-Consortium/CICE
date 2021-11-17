# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day
restart        gx3     4x2        debug,diag1
smoke          gbox80  1x1        box2001
smoke          gbox80  1x1        boxslotcyl
smoke	       gbox80  1x1        boxsymn
smoke	       gbox80  1x1        boxsyme
smoke	       gbox80  1x1        boxsymne

smoke          gx3     8x2        diag1,run5day,gridcd
restart        gx3     4x2        debug,diag1,gridcd
smoke          gbox80  1x1        box2001,gridcd
smoke          gbox80  1x1        boxslotcyl,gridcd
smoke	       gbox80  1x1        boxsymn,gridcd
smoke	       gbox80  1x1        boxsyme,gridcd
smoke	       gbox80  1x1        boxsymne,gridcd
