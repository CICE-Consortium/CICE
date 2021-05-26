# Test         Grid    PEs        Sets    BFB-compare
# some iobinary configurations fail due to bathymetry netcdf file requirement, remove them
# iobinary cannot work with JRA55 because netcdf is turned off
restart        gx3     8x4        gx3ncarbulk,debug,histall,iobinary,precision8
#restart        gx3    12x2        gx3ncarbulk,alt01,histall,iobinary
restart        gx3    16x2        gx3ncarbulk,alt02,histall,iobinary,precision8
#restart        gx3     4x2        gx3ncarbulk,alt03,histall,iobinary
restart        gx3     8x4        gx3ncarbulk,alt04,histall,iobinary,precision8
restart        gx3     4x4        gx3ncarbulk,alt05,histall,iobinary
restart        gx3    14x2        gx3ncarbulk,alt06,histall,iobinary,precision8
restart        gx3    32x1        gx3ncarbulk,bgcz,histall,iobinary,precision8
restart        gx3    16x2        gx3ncarbulk,bgcskl,histall,iobinary
restart        gx3    14x2        gx3ncarbulk,isotope,histall,iobinary,precision8
restart        gx3    16x2        gx3ncarbulk,fsd12,histall,iobinary

restart        gx3    32x1        debug,histall,ionetcdf
restart        gx3    15x2        alt01,histall,ionetcdf,precision8,cdf64
restart        gx3    15x2        alt02,histall,ionetcdf
restart        gx3    24x1        alt03,histall,ionetcdf,precision8
restart        gx3     8x4        alt04,histall,ionetcdf,cdf64
restart        gx3     8x4        alt05,histall,ionetcdf,precision8,cdf64
restart        gx3    16x2        alt06,histall,ionetcdf
restart        gx3    30x1        bgcz,histall,ionetcdf
restart        gx3    15x2        bgcskl,histall,ionetcdf,precision8
restart        gx3    31x1        isotope,histall,ionetcdf,cdf64
restart        gx3    14x2        fsd12,histall,ionetcdf,precision8

restart        gx3    16x2        debug,histall,iopio1,precision8,cdf64
restart        gx3    14x2        alt01,histall,iopio1,cdf64
restart        gx3    32x1        alt02,histall,iopio1,precision8
restart        gx3    24x1        alt03,histall,iopio1
restart        gx3     8x4        alt04,histall,iopio1,precision8,cdf64
restart        gx3     8x4        alt05,histall,iopio1,cdf64
restart        gx3    32x1        alt06,histall,iopio1,precision8
restart        gx3    16x2        bgcz,histall,iopio1,precision8
restart        gx3    30x1        bgcskl,histall,iopio1
restart        gx3     8x4        isotope,histall,iopio1,precision8,cdf64
restart        gx3    12x2        fsd12,histall,iopio1,cdf64

restart        gx3    16x2        debug,histall,iopio2
restart        gx3    14x2        alt01,histall,iopio2,precision8,cdf64
restart        gx3    32x1        alt02,histall,iopio2,cdf64
restart        gx3    24x1        alt03,histall,iopio2,precision8
restart        gx3     8x4        alt04,histall,iopio2
restart        gx3     8x4        alt05,histall,iopio2,precision8,cdf64
restart        gx3    16x2        alt06,histall,iopio2,cdf64
restart        gx3    16x2        bgcz,histall,iopio2,cdf64
restart        gx3    30x1        bgcskl,histall,iopio2,precision8
restart        gx3     8x4        isotope,histall,iopio2
restart        gx3    12x2        fsd12,histall,iopio2,precision8,cdf64

restart        gx3    16x2        debug,histall,iopio1p,precision8
restart        gx3    14x2        alt01,histall,iopio1p
restart        gx3    32x1        alt02,histall,iopio1p,precision8,cdf64
restart        gx3    24x1        alt03,histall,iopio1p,cdf64
restart        gx3     8x4        alt04,histall,iopio1p,precision8
restart        gx3     8x4        alt05,histall,iopio1p
restart        gx3     6x4        alt06,histall,iopio1p,precision8,cdf64
restart        gx3    16x2        bgcz,histall,iopio1p,precision8,cdf64
restart        gx3    30x1        bgcskl,histall,iopio1p,cdf64
restart        gx3     8x4        isotope,histall,iopio1p,precision8
restart        gx3    12x2        fsd12,histall,iopio1p

restart        gx3    16x2        debug,histall,iopio2p,cdf64
restart        gx3    14x2        alt01,histall,iopio2p,precision8
restart        gx3    32x1        alt02,histall,iopio2p
restart        gx3    24x1        alt03,histall,iopio2p,precision8,cdf64
restart        gx3     8x4        alt04,histall,iopio2p,cdf64
restart        gx3     8x4        alt05,histall,iopio2p,precision8
restart        gx3    24x1        alt06,histall,iopio2p
restart        gx3    16x2        bgcz,histall,iopio2p
restart        gx3    30x1        bgcskl,histall,iopio2p,precision8,cdf64
restart        gx3     8x4        isotope,histall,iopio2p,cdf64
restart        gx3    12x2        fsd12,histall,iopio2p,precision8

