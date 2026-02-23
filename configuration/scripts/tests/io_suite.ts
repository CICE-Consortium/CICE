# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3    32x1        run1year,histhrly,ionetcdf,iocdf2,short
smoke          gx3    32x1        run1year,histhrly,iopio2,iocdf5,short

# Need to check with cprnc manually, CICE tools always produce diffs between netcdf and pio
smoke          gx3     8x2        histall,ionetcdf,iocdf5,run5day
#smoke          gx3     8x2        histall,iopio2,iocdf5   smoke_gx3_8x2_histall_iocdf5_ionetcdf
smoke          gx3     8x2        histall,iopio2,iocdf5,run5day

# history restart tests
restart        gx3    15x2        gx3ncarbulk,fsd12,isotope,debug,histall10d,iobinary
restart        gx3    18x1        debug,fsd12,isotope,bgczm,histall10d,ionetcdf,iocdf5
restart        gx3    20x2        debug,fsd12,isotope,bgczm,histall10d,iopio1,iocdf5
restart        gx3    18x2        debug,fsd12,isotope,bgczm,histall10d,iopio2,iocdf2
restart        gx3    10x2        fsd12,isotope,bgczm,histall10d,ionetcdf,iocdf2
restart        gx3    40x1        fsd12,isotope,bgczm,histall10d,iopio1,iocdf1
restart        gx3    17x2        fsd12,isotope,bgczm,histall10d,iopio2,iocdf5

# some iobinary configurations fail due to bathymetry netcdf file requirement, remove them
# iobinary cannot work with JRA55 because netcdf is turned off
restart        gx3     8x4        gx3ncarbulk,debug,histall,iobinary,precision8
#restart        gx3    12x2        gx3ncarbulk,alt01,histall,iobinary
restart        gx3    16x2        gx3ncarbulk,alt02,histall,iobinary,precision8
#restart        gx3     4x2        gx3ncarbulk,alt03,histall,iobinary
restart        gx3     8x4        gx3ncarbulk,alt04,histall,iobinary,precision8
restart        gx3     4x4        gx3ncarbulk,alt05,histall,iobinary
restart        gx3    14x2        gx3ncarbulk,alt06,histall,iobinary,precision8
restart        gx3    14x2        gx3ncarbulk,alt07,histall,iobinary,precision8
#restart        gx3    32x1        gx3ncarbulk,bgcz,histall,iobinary,precision8
#restart        gx3    16x2        gx3ncarbulk,bgczm,histall,iobinary
restart        gx3    16x2        gx3ncarbulk,zaero,histall,iobinary
restart        gx3    14x2        gx3ncarbulk,isotope,histall,iobinary,precision8
restart        gx3    16x2        gx3ncarbulk,fsd12,histall,iobinary
restart        gx3     8x4        gx3ncarbulk,debug,histall,iobinary,precision8,histinst

restart        gx3    32x1        debug,histall,ionetcdf,iocdf1,precision8
restart        gx3    15x2        alt01,histall,ionetcdf,iocdf2,precision8
restart        gx3    15x2        alt02,histall,ionetcdf,iocdf5
restart        gx3    24x1        alt03,histall,ionetcdf,iohdf5,iohdf5opts
restart        gx3     8x4        alt04,histall,ionetcdf,iocdf1
restart        gx3     8x4        alt05,histall,ionetcdf,iocdf2
restart        gx3    16x2        alt06,histall,ionetcdf,iocdf5,precision8
restart        gx3    16x2        alt07,histall,ionetcdf,iohdf5,precision8
restart        gx3    30x1        bgczm,histall,ionetcdf,iocdf1
restart        gx3    15x2        bgcz,histall,ionetcdf,iocdf2,precision8
restart        gx3    31x1        isotope,histall,ionetcdf,iocdf5,precision8
restart        gx3    14x2        fsd12,histall,ionetcdf,iohdf5
restart        gx3    32x1        debug,histall,ionetcdf,iohdf5,histinst

restart        gx3    16x2x100x2x4  histall,iopio1,iopioopts
restart        gx3    16x2        debug,histall,iopio1,iocdf2
restart        gx3    14x2        alt01,histall,iopio1,iocdf5
restart        gx3    32x1        alt02,histall,iopio1,iohdf5
restart        gx3    24x1        alt03,histall,iopio1,iopnetcdf1,precision8
restart        gx3     8x4        alt04,histall,iopio1,iopnetcdf2,precision8
restart        gx3     8x4        alt05,histall,iopio1,iopnetcdf5,precision8
restart        gx3    32x1        alt06,histall,iopio1,iocdf1
restart        gx3    32x1        alt07,histall,iopio1,iocdf2,precision8
restart        gx3    16x2        bgczm,histall,iopio1,iocdf5,precision8
restart        gx3    30x1        bgcz,histall,iopio1,iohdf5,precision8
restart        gx3     8x4        isotope,histall,iopio1,iopnetcdf1
restart        gx3    12x2        fsd12,histall,iopio1,iopnetcdf2
restart        gx3    16x2        debug,histall,iopio1,iopnetcdf5,histinst

restart        gx3    16x2x100x2x4  debug,histall,iopio2,iopioopts
restart        gx3    16x2        debug,histall,iopio2,iopnetcdf1,precision8
restart        gx3    14x2        alt01,histall,iopio2,iopnetcdf2,precision8
restart        gx3    32x1        alt02,histall,iopio2,iopnetcdf5,precision8
restart        gx3    24x1        alt03,histall,iopio2,iocdf1
restart        gx3     8x4        alt04,histall,iopio2,iocdf2
restart        gx3     8x4        alt05,histall,iopio2,iocdf5
restart        gx3    16x2        alt06,histall,iopio2,iohdf5,iohdf5opts
restart        gx3    16x2        alt07,histall,iopio2,iopnetcdf1
restart        gx3    16x2        bgczm,histall,iopio2,iopnetcdf2
restart        gx3    30x1        bgcz,histall,iopio2,iopnetcdf5
restart        gx3     8x4        isotope,histall,iopio2,iohdf5,precision8
restart        gx3    12x2        fsd12,histall,iopio2,iocdf1,precision8
restart        gx3    16x2        debug,histall,iopio2,iocdf2,histinst,precision8

restart        gx3    12x2        cmip,ionetcdf,iocdf2
