# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day
restart        gx3     4x2        debug,diag1
restart2       gx1     16x2       debug,diag1
smoke          gbox12  1x1x12x12x1  boxchan
smoke          gbox80  1x1        box2001
smoke          gbox80  2x2        boxwallp5
smoke          gbox80  3x3        boxwall
smoke          gbox80  2x2        boxwallblock
smoke          gbox80  1x1        boxslotcyl
smoke          gbox80  2x4        boxnodyn
smoke	       gbox80  2x2        boxsymn
smoke	       gbox80  4x2        boxsyme
smoke	       gbox80  4x1        boxsymne
smoke	       gbox80  2x2        boxsymn,kmtislands
smoke	       gbox80  4x1        boxsyme,kmtislands
smoke	       gbox80  4x2        boxsymne,kmtislands
smoke	       gbox80  8x1        boxislandsn
smoke	       gbox80  4x2        boxislandse
smoke	       gbox80  2x4        boxislandsne
smoke          gx3     1x1x100x116x1   reprosum,run10day
smoke          gx3     1x1x25x29x16    reprosum,run10day,dwblockall  smoke_gx3_1x1x100x116x1_reprosum_run10day
smoke          gx3     1x1x5x4x580     reprosum,run10day,dwblockall  smoke_gx3_1x1x100x116x1_reprosum_run10day
smoke          gx3     1x1x5x4x580     reprosum,run10day             smoke_gx3_1x1x100x116x1_reprosum_run10day
smoke          gx1     32x1x16x16x32   reprosum,run10day
smoke          gx1     32x1x16x16x32   reprosum,run10day,cmplogrest,dwblockall   smoke_gx1_32x1x16x16x32_reprosum_run10day
smoke          gx1     32x1x16x12x40   reprosum,run10day,cmplogrest,dwblockall   smoke_gx1_32x1x16x16x32_reprosum_run10day
smoke          gx1     32x1x16x12x40   reprosum,run10day,cmplogrest              smoke_gx1_32x1x16x16x32_reprosum_run10day

smoke          gx3     8x2        diag1,run5day,gridcd
restart        gx3     4x2        debug,diag1,gridcd
restart2       gx1     16x2       debug,diag1,gridcd
smoke          gbox12  1x1x12x12x1  boxchan,gridcd
smoke          gbox80  1x1        box2001,gridcd
smoke          gbox80  2x2        boxwallp5,gridcd
smoke          gbox80  3x3        boxwall,gridcd
smoke          gbox80  2x2        boxwallblock,gridcd
smoke          gbox80  1x1        boxslotcyl,gridcd
smoke          gbox80  2x4        boxnodyn,gridcd
smoke	       gbox80  2x2        boxsymn,gridcd
smoke	       gbox80  4x2        boxsyme,gridcd
smoke	       gbox80  4x1        boxsymne,gridcd
smoke	       gbox80  2x2        boxsymn,kmtislands,gridcd
smoke	       gbox80  4x1        boxsyme,kmtislands,gridcd
smoke	       gbox80  4x2        boxsymne,kmtislands,gridcd
smoke	       gbox80  8x1        boxislandsn,gridcd
smoke	       gbox80  4x2        boxislandse,gridcd
smoke	       gbox80  2x4        boxislandsne,gridcd
smoke          gx3     1x1x100x116x1   reprosum,run10day,gridcd
smoke          gx3     1x1x25x29x16    reprosum,run10day,dwblockall,gridcd  smoke_gx3_1x1x100x116x1_gridcd_reprosum_run10day
smoke          gx3     1x1x5x4x580     reprosum,run10day,dwblockall,gridcd  smoke_gx3_1x1x100x116x1_gridcd_reprosum_run10day
smoke          gx3     1x1x5x4x580     reprosum,run10day,gridcd             smoke_gx3_1x1x100x116x1_gridcd_reprosum_run10day
smoke          gx1     32x1x16x16x32   reprosum,run10day,gridcd
smoke          gx1     32x1x16x16x32   reprosum,run10day,cmplogrest,dwblockall,gridcd   smoke_gx1_32x1x16x16x32_gridcd_reprosum_run10day
smoke          gx1     32x1x16x12x40   reprosum,run10day,cmplogrest,dwblockall,gridcd   smoke_gx1_32x1x16x16x32_gridcd_reprosum_run10day
smoke          gx1     32x1x16x12x40   reprosum,run10day,cmplogrest,gridcd              smoke_gx1_32x1x16x16x32_gridcd_reprosum_run10day

smoke          gx3     8x2        diag1,run5day,gridc
restart        gx3     4x2        debug,diag1,gridc
restart2       gx1     16x2       debug,diag1,gridc
smoke          gbox12  1x1x12x12x1  boxchan,gridc
smoke          gbox80  1x1        box2001,gridc
smoke          gbox80  2x2        boxwallp5,gridc
smoke          gbox80  3x3        boxwall,gridc
smoke          gbox80  2x2        boxwallblock,gridc
smoke          gbox80  1x1        boxslotcyl,gridc
smoke          gbox80  2x4        boxnodyn,gridc
smoke	       gbox80  2x2        boxsymn,gridc
smoke	       gbox80  4x2        boxsyme,gridc
smoke	       gbox80  4x1        boxsymne,gridc
smoke	       gbox80  2x2        boxsymn,kmtislands,gridc
smoke	       gbox80  4x1        boxsyme,kmtislands,gridc
smoke	       gbox80  4x2        boxsymne,kmtislands,gridc
smoke	       gbox80  8x1        boxislandsn,gridc
smoke	       gbox80  4x2        boxislandse,gridc
smoke	       gbox80  2x4        boxislandsne,gridc
smoke          gx3     1x1x100x116x1   reprosum,run10day,gridc
smoke          gx3     1x1x25x29x16    reprosum,run10day,dwblockall,gridc  smoke_gx3_1x1x100x116x1_gridc_reprosum_run10day
smoke          gx3     1x1x5x4x580     reprosum,run10day,dwblockall,gridc  smoke_gx3_1x1x100x116x1_gridc_reprosum_run10day
smoke          gx3     1x1x5x4x580     reprosum,run10day,gridc             smoke_gx3_1x1x100x116x1_gridc_reprosum_run10day
smoke          gx1     32x1x16x16x32   reprosum,run10day,gridc
smoke          gx1     32x1x16x16x32   reprosum,run10day,cmplogrest,dwblockall,gridc   smoke_gx1_32x1x16x16x32_gridc_reprosum_run10day
smoke          gx1     32x1x16x12x40   reprosum,run10day,cmplogrest,dwblockall,gridc   smoke_gx1_32x1x16x16x32_gridc_reprosum_run10day
smoke          gx1     32x1x16x12x40   reprosum,run10day,cmplogrest,gridc              smoke_gx1_32x1x16x16x32_gridc_reprosum_run10day
