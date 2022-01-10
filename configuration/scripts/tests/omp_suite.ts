# Test         Grid    PEs        Sets    BFB-compare

smoke          gx3     8x4        diag1,reprosum,run10day
smoke          gx3     6x2        alt01,reprosum,run10day
smoke          gx3     8x2        alt02,reprosum,run10day
smoke          gx3    12x2        alt03,droundrobin,reprosum,run10day
smoke          gx3     4x4        alt04,reprosum,run10day
smoke          gx3     4x4        alt05,reprosum,run10day
smoke          gx3     8x2        alt06,reprosum,run10day
smoke          gx3     8x2        bgcz,reprosum,run10day
smoke          gx1    15x2        seabedprob,reprosum,run10day
smoke          gx3    14x2        fsd12,reprosum,run10day
smoke          gx3    11x2        isotope,reprosum,run10day
smoke          gx3     8x4        snwitdrdg,snwgrain,icdefault,reprosum,run10day
smoke          gx3     6x4        dynpicard,reprosum,run10day
smoke          gx3     8x3        zsal,reprosum,run10day

smoke        gbox128   8x2        reprosum,run10day
smoke        gbox128  12x2        boxnodyn,reprosum,run10day
smoke        gbox128   9x2        boxadv,reprosum,run10day
smoke        gbox128  14x2        boxrestore,reprosum,run10day
smoke        gbox80    4x5        box2001,reprosum,run10day
smoke        gbox80   11x3        boxslotcyl,reprosum,run10day

smoke          gx3     4x2        diag1,reprosum,run10day,cmplogrest              smoke_gx3_8x4_diag1_reprosum_run10day
smoke          gx3     4x1        diag1,reprosum,run10day,cmplogrest,thread       smoke_gx3_8x4_diag1_reprosum_run10day
smoke          gx3     8x1        alt01,reprosum,run10day,cmplogrest,thread       smoke_gx3_6x2_alt01_reprosum_run10day
smoke          gx3    16x1        alt02,reprosum,run10day,cmplogrest,thread       smoke_gx3_8x2_alt02_reprosum_run10day
smoke          gx3    24x1        alt03,reprosum,run10day,cmplogrest,thread       smoke_gx3_12x2_alt03_droundrobin_reprosum_run10day
smoke          gx3    24x1        alt04,reprosum,run10day,cmplogrest,thread       smoke_gx3_4x4_alt04_reprosum_run10day
smoke          gx3    14x1        alt05,reprosum,run10day,cmplogrest,thread       smoke_gx3_4x4_alt05_reprosum_run10day
smoke          gx3    24x1        alt06,reprosum,run10day,cmplogrest,thread       smoke_gx3_8x2_alt06_reprosum_run10day
smoke          gx3    12x1        bgcz,reprosum,run10day,cmplogrest,thread        smoke_gx3_8x2_bgcz_reprosum_run10day
smoke          gx1    28x1        seabedprob,reprosum,run10day,cmplogrest,thread  smoke_gx1_15x2_reprosum_run10day_seabedprob
smoke          gx3    30x1        fsd12,reprosum,run10day,cmplogrest,thread       smoke_gx3_14x2_fsd12_reprosum_run10day
smoke          gx3    16x1        isotope,reprosum,run10day,cmplogrest,thread     smoke_gx3_11x2_isotope_reprosum_run10day
smoke          gx3    28x1        snwitdrdg,snwgrain,icdefault,reprosum,run10day,cmplogrest,thread  smoke_gx3_8x4_icdefault_reprosum_run10day_snwitdrdg_snwgrain
smoke          gx3    18x1        dynpicard,reprosum,run10day,cmplogrest,thread   smoke_gx3_6x4_dynpicard_reprosum_run10day
smoke          gx3    20x1        zsal,reprosum,run10day,cmplogrest,thread        smoke_gx3_8x3_reprosum_run10day_zsal

smoke        gbox128  20x1        reprosum,run10day,cmplogrest,thread             smoke_gbox128_8x2_reprosum_run10day
smoke        gbox128  16x1        boxnodyn,reprosum,run10day,cmplogrest,thread    smoke_gbox128_12x2_boxnodyn_reprosum_run10day
smoke        gbox128  14x1        boxadv,reprosum,run10day,cmplogrest,thread      smoke_gbox128_9x2_boxadv_reprosum_run10day
smoke        gbox128  24x1        boxrestore,reprosum,run10day,cmplogrest,thread  smoke_gbox128_14x2_boxrestore_reprosum_run10day
smoke        gbox80   19x1        box2001,reprosum,run10day,cmplogrest,thread     smoke_gbox80_4x5_box2001_reprosum_run10day
smoke        gbox80    8x4        boxslotcyl,reprosum,run10day,cmplogrest,thread  smoke_gbox80_11x3_boxslotcyl_reprosum_run10day
