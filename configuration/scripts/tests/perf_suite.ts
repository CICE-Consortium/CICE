# Test         Grid    PEs        Sets    BFB-compare
smoke        gx1      1x1x320x384x1  run2day,droundrobin
smoke        gx1     64x1x16x16x8    run2day,droundrobin,thread
sleep 180
#
smoke        gx1      1x1x320x384x1  run2day,droundrobin
smoke        gx1      1x1x160x192x4  run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1      1x1x80x96x16   run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1      1x1x40x48x64   run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1      1x1x20x24x256  run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
#
smoke        gx1      1x1x16x16x480  run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1      2x1x16x16x240  run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1      4x1x16x16x120  run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1      8x1x16x16x60   run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1     16x1x16x16x30   run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1     32x1x16x16x15   run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1     64x1x16x16x8    run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1    128x1x16x16x4    run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
#
smoke        gx1     64x1x16x16x8    run2day,droundrobin  smoke_gx1_1x1x320x384x1_droundrobin_run2day
smoke        gx1     64x1x16x16x8    run2day,droundrobin,thread
smoke        gx1     32x2x16x16x16   run2day,droundrobin  smoke_gx1_64x1x16x16x8_droundrobin_run2day_thread
smoke        gx1     16x4x16x16x32   run2day,droundrobin  smoke_gx1_64x1x16x16x8_droundrobin_run2day_thread
smoke        gx1      8x8x16x16x64   run2day,droundrobin  smoke_gx1_64x1x16x16x8_droundrobin_run2day_thread
smoke        gx1     4x16x16x16x128  run2day,droundrobin  smoke_gx1_64x1x16x16x8_droundrobin_run2day_thread
smoke        gx1     32x2x16x16x16   run2day,droundrobin,ompscheds   smoke_gx1_64x1x16x16x8_droundrobin_run2day_thread
smoke        gx1     32x2x16x16x16   run2day,droundrobin,ompschedd1  smoke_gx1_64x1x16x16x8_droundrobin_run2day_thread
smoke        gx1     32x2x16x16x16   run2day,droundrobin,ompscheds1  smoke_gx1_64x1x16x16x8_droundrobin_run2day_thread
#
