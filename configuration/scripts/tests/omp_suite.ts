# Test         Grid    PEs        Sets    BFB-compare

smoke          gx3     8x4        diag1,reprosum,run10day
smoke          gx3     6x2        alt01,reprosum,run10day
smoke          gx3     8x2        alt02,reprosum,run10day
smoke          gx3    12x2        alt03,droundrobin,reprosum,run10day
smoke          gx3     4x4        alt04,reprosum,run10day
smoke          gx3     4x4        alt05,reprosum,run10day
smoke          gx3     8x2        alt06,reprosum,run10day
smoke          gx3     8x2        bgcz,reprosum,run10day
smoke          gx1    15x2        reprosum,run10day
smoke          gx1    15x2        seabedprob,reprosum,run10day
smoke          gx3    14x2        fsd12,reprosum,run10day
smoke          gx3    11x2        isotope,reprosum,run10day
smoke          gx3     8x4        snwitdrdg,snwgrain,icdefault,reprosum,run10day
smoke          gx3     6x4        dynpicard,reprosum,run10day
smoke          gx3     8x3        zsal,reprosum,run10day
smoke          gx3     1x1x100x116x1   reprosum,run10day,thread

smoke        gbox128   8x2        reprosum,run10day
smoke        gbox128  12x2        boxnodyn,reprosum,run10day
smoke        gbox128   9x2        boxadv,reprosum,run10day
smoke        gbox128  14x2        boxrestore,reprosum,run10day
smoke        gbox80    4x5        box2001,reprosum,run10day
smoke        gbox80   11x3        boxslotcyl,reprosum,run10day

smoke          gx3     4x2        diag1,reprosum,run10day,cmplogrest              smoke_gx3_8x4_diag1_reprosum_run10day
smoke          gx3     4x1        diag1,reprosum,run10day,cmplogrest,thread       smoke_gx3_8x4_diag1_reprosum_run10day
smoke          gx3     8x1        alt01,reprosum,run10day,cmplogrest,thread       smoke_gx3_6x2_alt01_reprosum_run10day
smoke          gx3     8x1        alt02,reprosum,run10day,cmplogrest,thread       smoke_gx3_8x2_alt02_reprosum_run10day
smoke          gx3     8x1        alt03,reprosum,run10day,cmplogrest,thread       smoke_gx3_12x2_alt03_droundrobin_reprosum_run10day
smoke          gx3     8x1        alt04,reprosum,run10day,cmplogrest,thread       smoke_gx3_4x4_alt04_reprosum_run10day
smoke          gx3     8x1        alt05,reprosum,run10day,cmplogrest,thread       smoke_gx3_4x4_alt05_reprosum_run10day
smoke          gx3     8x1        alt06,reprosum,run10day,cmplogrest,thread       smoke_gx3_8x2_alt06_reprosum_run10day
smoke          gx3     8x1        bgcz,reprosum,run10day,cmplogrest,thread        smoke_gx3_8x2_bgcz_reprosum_run10day
smoke          gx1    18x1        reprosum,run10day,cmplogrest,thread             smoke_gx1_15x2_reprosum_run10day
smoke          gx1    18x1        seabedprob,reprosum,run10day,cmplogrest,thread  smoke_gx1_15x2_reprosum_run10day_seabedprob
smoke          gx3     8x1        fsd12,reprosum,run10day,cmplogrest,thread       smoke_gx3_14x2_fsd12_reprosum_run10day
smoke          gx3     8x1        isotope,reprosum,run10day,cmplogrest,thread     smoke_gx3_11x2_isotope_reprosum_run10day
smoke          gx3     8x1        snwitdrdg,snwgrain,icdefault,reprosum,run10day,cmplogrest,thread  smoke_gx3_8x4_icdefault_reprosum_run10day_snwitdrdg_snwgrain
smoke          gx3     8x1        dynpicard,reprosum,run10day,cmplogrest,thread   smoke_gx3_6x4_dynpicard_reprosum_run10day
smoke          gx3     8x1        zsal,reprosum,run10day,cmplogrest,thread        smoke_gx3_8x3_reprosum_run10day_zsal
smoke          gx3     4x2x25x29x4     reprosum,run10day                          smoke_gx3_1x1x100x116x1_reprosum_run10day_thread
smoke          gx3     8x4x5x4x80      reprosum,run10day                          smoke_gx3_1x1x100x116x1_reprosum_run10day_thread

smoke        gbox128   8x1        reprosum,run10day,cmplogrest,thread             smoke_gbox128_8x2_reprosum_run10day
smoke        gbox128   8x1        boxnodyn,reprosum,run10day,cmplogrest,thread    smoke_gbox128_12x2_boxnodyn_reprosum_run10day
smoke        gbox128   8x1        boxadv,reprosum,run10day,cmplogrest,thread      smoke_gbox128_9x2_boxadv_reprosum_run10day
smoke        gbox128   8x1        boxrestore,reprosum,run10day,cmplogrest,thread  smoke_gbox128_14x2_boxrestore_reprosum_run10day
smoke        gbox80    8x1        box2001,reprosum,run10day,cmplogrest,thread     smoke_gbox80_4x5_box2001_reprosum_run10day
smoke        gbox80    8x1        boxslotcyl,reprosum,run10day,cmplogrest,thread  smoke_gbox80_11x3_boxslotcyl_reprosum_run10day

#gridC

smoke          gx3     8x4        diag1,reprosum,run10day,gridc
smoke          gx3     6x2        alt01,reprosum,run10day,gridc
smoke          gx3     8x2        alt02,reprosum,run10day,gridc
#smoke          gx3    12x2        alt03,droundrobin,reprosum,run10day,gridc
smoke          gx3     4x4        alt04,reprosum,run10day,gridc
smoke          gx3     4x4        alt05,reprosum,run10day,gridc
smoke          gx3     8x2        alt06,reprosum,run10day,gridc
smoke          gx3     8x2        bgcz,reprosum,run10day,gridc
smoke          gx1    15x2        reprosum,run10day,gridc
smoke          gx1    15x2        seabedprob,reprosum,run10day,gridc
smoke          gx3    14x2        fsd12,reprosum,run10day,gridc
smoke          gx3    11x2        isotope,reprosum,run10day,gridc
smoke          gx3     8x4        snwitdrdg,snwgrain,icdefault,reprosum,run10day,gridc
#smoke          gx3     6x4        dynpicard,reprosum,run10day,gridc
smoke          gx3     8x3        zsal,reprosum,run10day,gridc
smoke          gx3     1x1x100x116x1   reprosum,run10day,gridc,thread

smoke        gbox128   8x2        reprosum,run10day,gridc
smoke        gbox128  12x2        boxnodyn,reprosum,run10day,gridc
#smoke        gbox128   9x2        boxadv,reprosum,run10day,gridc
smoke        gbox128  14x2        boxrestore,reprosum,run10day,gridc
smoke        gbox80    4x5        box2001,reprosum,run10day,gridc
smoke        gbox80   11x3        boxslotcyl,reprosum,run10day,gridc

smoke          gx3     4x2        diag1,reprosum,run10day,cmplogrest,gridc              smoke_gx3_8x4_gridc_diag1_reprosum_run10day
smoke          gx3     4x1        diag1,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_8x4_gridc_diag1_reprosum_run10day
smoke          gx3     8x1        alt01,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_6x2_alt01_gridc_reprosum_run10day
smoke          gx3     8x1        alt02,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_8x2_alt02_gridc_reprosum_run10day
#smoke          gx3     8x1        alt03,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_12x2_alt03_droundrobin_gridc_reprosum_run10day
smoke          gx3     8x1        alt04,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_4x4_alt04_gridc_reprosum_run10day
smoke          gx3     8x1        alt05,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_4x4_alt05_gridc_reprosum_run10day
smoke          gx3     8x1        alt06,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_8x2_alt06_gridc_reprosum_run10day
smoke          gx3     8x1        bgcz,reprosum,run10day,cmplogrest,thread,gridc        smoke_gx3_8x2_bgcz_gridc_reprosum_run10day
smoke          gx1    18x1        reprosum,run10day,cmplogrest,thread,gridc             smoke_gx1_15x2_gridc_reprosum_run10day
smoke          gx1    18x1        seabedprob,reprosum,run10day,cmplogrest,thread,gridc  smoke_gx1_15x2_gridc_reprosum_run10day_seabedprob
smoke          gx3     8x1        fsd12,reprosum,run10day,cmplogrest,thread,gridc       smoke_gx3_14x2_fsd12_gridc_reprosum_run10day
smoke          gx3     8x1        isotope,reprosum,run10day,cmplogrest,thread,gridc     smoke_gx3_11x2_gridc_isotope_reprosum_run10day
smoke          gx3     8x1        snwitdrdg,snwgrain,icdefault,reprosum,run10day,cmplogrest,thread,gridc  smoke_gx3_8x4_gridc_icdefault_reprosum_run10day_snwitdrdg_snwgrain
#smoke          gx3     8x1        dynpicard,reprosum,run10day,cmplogrest,thread,gridc   smoke_gx3_6x4_dynpicard_gridc_reprosum_run10day
smoke          gx3     8x1        zsal,reprosum,run10day,cmplogrest,thread,gridc        smoke_gx3_8x3_gridc_reprosum_run10day_zsal
smoke          gx3     4x2x25x29x4     reprosum,run10day,gridc                          smoke_gx3_1x1x100x116x1_gridc_reprosum_run10day_thread
smoke          gx3     8x4x5x4x80      reprosum,run10day,gridc                          smoke_gx3_1x1x100x116x1_gridc_reprosum_run10day_thread

smoke        gbox128   8x1        reprosum,run10day,cmplogrest,thread,gridc             smoke_gbox128_8x2_gridc_reprosum_run10day
smoke        gbox128   8x1        boxnodyn,reprosum,run10day,cmplogrest,thread,gridc    smoke_gbox128_12x2_boxnodyn_gridc_reprosum_run10day
#smoke        gbox128   8x1        boxadv,reprosum,run10day,cmplogrest,thread,gridc      smoke_gbox128_9x2_boxadv_gridc_reprosum_run10day
smoke        gbox128   8x1        boxrestore,reprosum,run10day,cmplogrest,thread,gridc  smoke_gbox128_14x2_boxrestore_gridc_reprosum_run10day
smoke        gbox80    8x1        box2001,reprosum,run10day,cmplogrest,thread,gridc     smoke_gbox80_4x5_box2001_gridc_reprosum_run10day
smoke        gbox80    8x1        boxslotcyl,reprosum,run10day,cmplogrest,thread,gridc  smoke_gbox80_11x3_boxslotcyl_gridc_reprosum_run10day

#gridCD

smoke          gx3     8x4        diag1,reprosum,run10day,gridcd
smoke          gx3     6x2        alt01,reprosum,run10day,gridcd
smoke          gx3     8x2        alt02,reprosum,run10day,gridcd
#smoke          gx3    12x2        alt03,droundrobin,reprosum,run10day,gridcd
smoke          gx3     4x4        alt04,reprosum,run10day,gridcd
smoke          gx3     4x4        alt05,reprosum,run10day,gridcd
smoke          gx3     8x2        alt06,reprosum,run10day,gridcd
smoke          gx3     8x2        bgcz,reprosum,run10day,gridcd
smoke          gx1    15x2        reprosum,run10day,gridcd
smoke          gx1    15x2        seabedprob,reprosum,run10day,gridcd
smoke          gx3    14x2        fsd12,reprosum,run10day,gridcd
smoke          gx3    11x2        isotope,reprosum,run10day,gridcd
smoke          gx3     8x4        snwitdrdg,snwgrain,icdefault,reprosum,run10day,gridcd
#smoke          gx3     6x4        dynpicard,reprosum,run10day,gridcd
smoke          gx3     8x3        zsal,reprosum,run10day,gridcd
smoke          gx3     1x1x100x116x1   reprosum,run10day,gridcd,thread

smoke        gbox128   8x2        reprosum,run10day,gridcd
smoke        gbox128  12x2        boxnodyn,reprosum,run10day,gridcd
#smoke        gbox128   9x2        boxadv,reprosum,run10day,gridcd
smoke        gbox128  14x2        boxrestore,reprosum,run10day,gridcd
smoke        gbox80    4x5        box2001,reprosum,run10day,gridcd
smoke        gbox80   11x3        boxslotcyl,reprosum,run10day,gridcd

smoke          gx3     4x2        diag1,reprosum,run10day,cmplogrest,gridcd              smoke_gx3_8x4_gridcd_diag1_reprosum_run10day
smoke          gx3     4x1        diag1,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_8x4_gridcd_diag1_reprosum_run10day
smoke          gx3     8x1        alt01,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_6x2_alt01_gridcd_reprosum_run10day
smoke          gx3     8x1        alt02,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_8x2_alt02_gridcd_reprosum_run10day
#smoke          gx3     8x1        alt03,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_12x2_alt03_droundrobin_gridcd_reprosum_run10day
smoke          gx3     8x1        alt04,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_4x4_alt04_gridcd_reprosum_run10day
smoke          gx3     8x1        alt05,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_4x4_alt05_gridcd_reprosum_run10day
smoke          gx3     8x1        alt06,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_8x2_alt06_gridcd_reprosum_run10day
smoke          gx3     8x1        bgcz,reprosum,run10day,cmplogrest,thread,gridcd        smoke_gx3_8x2_bgcz_gridcd_reprosum_run10day
smoke          gx1    18x1        reprosum,run10day,cmplogrest,thread,gridcd             smoke_gx1_15x2_gridcd_reprosum_run10day
smoke          gx1    18x1        seabedprob,reprosum,run10day,cmplogrest,thread,gridcd  smoke_gx1_15x2_gridcd_reprosum_run10day_seabedprob
smoke          gx3     8x1        fsd12,reprosum,run10day,cmplogrest,thread,gridcd       smoke_gx3_14x2_fsd12_gridcd_reprosum_run10day
smoke          gx3     8x1        isotope,reprosum,run10day,cmplogrest,thread,gridcd     smoke_gx3_11x2_gridcd_isotope_reprosum_run10day
smoke          gx3     8x1        snwitdrdg,snwgrain,icdefault,reprosum,run10day,cmplogrest,thread,gridcd  smoke_gx3_8x4_gridcd_icdefault_reprosum_run10day_snwitdrdg_snwgrain
#smoke          gx3     8x1        dynpicard,reprosum,run10day,cmplogrest,thread,gridcd   smoke_gx3_6x4_dynpicard_gridcd_reprosum_run10day
smoke          gx3     8x1        zsal,reprosum,run10day,cmplogrest,thread,gridcd        smoke_gx3_8x3_gridcd_reprosum_run10day_zsal
smoke          gx3     4x2x25x29x4     reprosum,run10day,gridcd                          smoke_gx3_1x1x100x116x1_gridcd_reprosum_run10day_thread
smoke          gx3     8x4x5x4x80      reprosum,run10day,gridcd                          smoke_gx3_1x1x100x116x1_gridcd_reprosum_run10day_thread

smoke        gbox128   8x1        reprosum,run10day,cmplogrest,thread,gridcd             smoke_gbox128_8x2_gridcd_reprosum_run10day
smoke        gbox128   8x1        boxnodyn,reprosum,run10day,cmplogrest,thread,gridcd    smoke_gbox128_12x2_boxnodyn_gridcd_reprosum_run10day
#smoke        gbox128   8x1        boxadv,reprosum,run10day,cmplogrest,thread,gridcd      smoke_gbox128_9x2_boxadv_gridcd_reprosum_run10day
smoke        gbox128   8x1        boxrestore,reprosum,run10day,cmplogrest,thread,gridcd  smoke_gbox128_14x2_boxrestore_gridcd_reprosum_run10day
smoke        gbox80    8x1        box2001,reprosum,run10day,cmplogrest,thread,gridcd     smoke_gbox80_4x5_box2001_gridcd_reprosum_run10day
smoke        gbox80    8x1        boxslotcyl,reprosum,run10day,cmplogrest,thread,gridcd  smoke_gbox80_11x3_boxslotcyl_gridcd_reprosum_run10day

