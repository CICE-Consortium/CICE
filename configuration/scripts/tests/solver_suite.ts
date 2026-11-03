# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     1x8        diag1,run5day
smoke          gx3     1x8        diag1,run5day,evp1d       smoke_gx3_1x8_diag1_run5day
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day,evp1d smoke_gbox80_8x1_boxforcee_boxopen_kmtislands_run1day
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day,boxridge
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day,boxridge,evp1d smoke_gbox80_8x1_boxforcee_boxopen_boxridge_kmtislands_run1day
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day,nghost2
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day,nghost2,evp1d smoke_gbox80_8x1_boxforcee_boxopen_kmtislands_nghost2_run1day
smoke          gx3     4x4        histall,run2day
smoke          gx3     4x4        histall,run2day,evp1d     smoke_gx3_4x4_histall_run2day
smoke          gx3     4x1        dynpicard
