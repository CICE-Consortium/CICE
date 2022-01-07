# Test         Grid    PEs        Sets    BFB-compare

smoke          gx3     8x4        diag1
smoke          gx3     6x2        alt01
smoke          gx3     8x2        alt02
smoke          gx3     12x2        alt03
smoke          gx3     4x4        alt04
smoke          gx3     4x4        alt05
smoke          gx3     8x2        alt06
smoke          gx3     8x2        bgcz
smoke          gx1    15x2        seabedprob
smoke          gx3    14x2        fsd12
smoke          gx3    11x2        isotope
smoke          gx3     8x4        snwitdrdg,snwgrain,icdefault
smoke          gx3     6x4        dynpicard
smoke          gx3     8x3        zsal

smoke        gbox128   8x2        short
smoke        gbox128  12x2        boxnodyn,short
smoke        gbox128   9x2        boxadv,short
smoke        gbox128  14x2        boxrestore,short
smoke        gbox80    4x5        box2001
smoke        gbox80   11x3        boxslotcyl

smoke          gx3     4x2        diag1              smoke_gx3_8x2_diag1
smoke          gx3     4x1        diag1,thread       smoke_gx3_8x2_diag1
smoke          gx3     8x1        alt01,thread       smoke_gx3_6x2_alt01
smoke          gx3    16x1        alt02,thread       smoke_gx3_8x2_alt02
smoke          gx3    24x1        alt03,thread       smoke_gx3_12x2_alt03
smoke          gx3    24x1        alt04,thread       smoke_gx3_4x4_alt04
smoke          gx3    14x1        alt05,thread       smoke_gx3_4x4_alt05
smoke          gx3    24x1        alt06,thread       smoke_gx3_8x2_alt06
smoke          gx3    12x1        bgcz,thread        smoke_gx3_8x2_bgcz
smoke          gx1    28x1        seabedprob,thread  smoke_gx1_15x2_seabedprob
smoke          gx3    30x1        fsd12,thread       smoke_gx3_14x2_fsd12
smoke          gx3    16x1        isotope,thread     smoke_gx3_11x2_isotope
smoke          gx3    28x1        snwitdrdg,snwgrain,icdefault,thread  smoke_gx3_8x4_icdefault_snwitdrdg_snwgrain
smoke          gx3    18x1        dynpicard,thread   smoke_gx3_6x4_dynpicard
smoke          gx3    20x1        zsal,thread        smoke_gx3_8x3_zsal

smoke        gbox128  20x1        short,thread            smoke_gbox128_8x2_short
smoke        gbox128  16x1        boxnodyn,short,thread   smoke_gbox128_12x2_boxnodyn_short
smoke        gbox128  14x1        boxadv,short,thread     smoke_gbox128_9x2_boxadv_short
smoke        gbox128  24x1        boxrestore,short,thread smoke_gbox128_14x2_boxrestore_short
smoke        gbox80   19x1        box2001,thread          smoke_gbox80_4x5_box2001
smoke        gbox80    8x4        boxslotcyl,thread       smoke_gbox80_11x3_boxslotcyl


