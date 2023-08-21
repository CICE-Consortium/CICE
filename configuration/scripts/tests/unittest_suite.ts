# Test         Grid    PEs          Sets         BFB-compare
smoke          gx3     8x2           diag1,run5day
smoke          gx3     4x2x25x29x4   debug,run2day,dslenderX2
unittest       gx3     1x1           helloworld
unittest       gx3     1x1           calchk,short
unittest       gx3     4x1x25x29x4   sumchk
unittest       gx3     1x1x25x29x16  sumchk
unittest       tx1     8x1           sumchk
unittest       tx1     8x1           sumchk,tripolet
unittest       gx3     4x1           bcstchk
unittest       gx3     1x1           bcstchk
unittest       gx3     8x2           gridavgchk,dwblockall
unittest       gx3     12x1          gridavgchk
unittest       gx1     28x1          gridavgchk,dwblockall
unittest       gx1     16x2          gridavgchk
unittest       gbox128 8x2           gridavgchk
unittest       gbox80  1x1x10x10x80  halochk,cyclic,debug
unittest       gbox80  1x1x24x23x16  halochk
unittest       gbox80  1x1x23x24x16  halochk,cyclic
unittest       gbox80  1x1x23x23x16  halochk,open
unittest       tx1     1x1x90x60x16  halochk,dwblockall
unittest       tx1     1x1x90x60x16  halochk,dwblockall,tripolet
unittest       tx1     1x1x95x65x16  halochk,dwblockall
unittest       tx1     1x1x95x65x16  halochk,dwblockall,tripolet
unittest       gx3     4x2           halochk,dwblockall,debug
unittest       gx3     8x2x16x12x10  halochk,cyclic,dwblockall
unittest       gx3     17x1x16x12x10 halochk,open,dwblockall
unittest       tx1     4x2           halochk,dwblockall
unittest       tx1     4x2           halochk,dwblockall,tripolet
unittest       tx1     4x2x65x45x10  halochk,dwblockall
unittest       tx1     4x2x57x43x12  halochk,dwblockall,tripolet
unittest       gx3     1x1           optargs
unittest       gx3     1x1           opticep
unittest       gx3     4x2x25x29x4   debug,run2day,dslenderX2,opticep,cmplog  smoke_gx3_4x2x25x29x4_debug_dslenderX2_run2day
unittest       gx3     8x2           diag1,run5day,opticep,cmplog             smoke_gx3_8x2_diag1_run5day
