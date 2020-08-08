# Test         Grid    PEs        Sets    BFB-compare
# some iobinary configurations fail due to bathymetry netcdf file requirement, remove them
restart        gx3     8x4        debug,histall,iobinary,precision8
#restart        gx3    12x2        alt01,histall,iobinary
restart        gx3    16x2        alt02,histall,iobinary,precision8
#restart        gx3     4x2        alt03,histall,iobinary
restart        gx3     8x4        alt04,histall,iobinary,precision8
restart        gx3     4x4        alt05,histall,iobinary
restart        gx3    32x1        bgcz,histall,iobinary,precision8
restart        gx3    16x2        bgcskl,histall,iobinary
restart        gx3    14x2        isotope,histall,iobinary,precision8
restart        gx3    16x2        fsd12,histall,iobinary

restart        gx3    32x1        debug,histall,ionetcdf
restart        gx3    15x2        alt01,histall,ionetcdf,precision8,cdf64
restart        gx3    15x2        alt02,histall,ionetcdf
restart        gx3    24x1        alt03,histall,ionetcdf,precision8
restart        gx3     8x4        alt04,histall,ionetcdf,cdf64
restart        gx3     8x4        alt05,histall,ionetcdf,precision8,cdf64
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
restart        gx3    16x2        bgcz,histall,iopio2p
restart        gx3    30x1        bgcskl,histall,iopio2p,precision8,cdf64
restart        gx3     8x4        isotope,histall,iopio2p,cdf64
restart        gx3    12x2        fsd12,histall,iopio2p,precision8

