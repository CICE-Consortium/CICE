# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day
smoke          gx3     1x1        debug,diag1,run2day
smoke          gx3     1x4        debug,diag1,run2day
smoke          gx3     4x1        debug,diag1,run5day
restart        gx3     8x2        debug
restart        gx3     8x2        debug,gx3nc
smoke          gx3     8x2        diag24,run1year,long
smoke          gx3     7x2        diag1,bigdiag,run1day,diagpt1
decomp         gx3     4x2x25x29x5  none
smoke          gx3     4x2        diag1,run5day             smoke_gx3_8x2_diag1_run5day
smoke          gx3     4x1        diag1,run5day,thread      smoke_gx3_8x2_diag1_run5day
# EVP 1d vs 2d.  Nothing in the suite had ever compared the two solvers against
# each other, and that gap is how stress_eap stayed wrong until c6df487 and how
# the MOM-supergrid DminTarea divergence survived until 03bc112: a restart test
# compares a run to itself, so a systematic error cancels.  Each pair below is
# the same grid, PE layout and length, differing only in evp_algorithm.
#
# The debug pair is the definitive one.  At -O0 neither path is vectorised or
# reassociated, so a difference is algorithmic.  The optimised pair can differ
# by an ULP purely because the two code paths vectorise differently, so read it
# second.
#
# NB: these are meaningless in a sweep that forces evp_algorithm globally --
# both halves would then run the 1d solver and compare it with itself.
smoke          gx3     1x4        debug,diag1,run2day,evp1d smoke_gx3_1x4_debug_diag1_run2day
smoke          gx3     1x8        diag1,run5day
smoke          gx3     1x8        diag1,run5day,evp1d       smoke_gx3_1x8_diag1_run5day
#
# A box-grid pair.  The pairs above are all gx3; this adds a second grid, and
# gbox80 is cheap.  Verified to compare clean on main before being added.
#
# NB: the same pair with bclinearextrap is NOT here, deliberately.  On main it
# fails twice over, and neither is ours: the evp1d half aborts under bounds
# checking -- convert_2d_1d_init builds the north-neighbour indices Isw/Isse
# as i(+1) + (j-0)*nx with no guard for j = ny, so a top-row cell indexes one
# row past the extended grid and convert_2d_1d_dyn reads G_uvel(i,ny+1) at
# ice_dyn_evp1d.F90:947 -- and the two solvers then disagree, plausibly
# because that out-of-bounds read lands in the halo velocities.  Plain
# boundaries put no active T cell in the top row, so neither shows up here.
# Add that pair once the index arithmetic is fixed upstream.
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day
smoke          gbox80  8x1        boxopen,kmtislands,boxforcee,run1day,evp1d smoke_gbox80_8x1_boxopen_kmtislands_boxforcee_run1day
#
# And one pair that writes output while running the 1d solver.  evp1d scatters
# the twelve stresses back only when a file is due, so this reaches a path the
# two pairs above never take.
#
# Be clear about what it does and does not prove.  comparebfb compares iced*
# restart files only -- no history file is bit-compared anywhere in CICE, and
# they carry creation timestamps so they never could be.  So the restart half
# of the gate is verified: stale stresses would land in iced* and this would
# fail.  The history half runs but its output is not checked, and a fault
# confined to it would be invisible to the whole suite.  That is a limit of the
# machinery, not of this test; keep the gate condition simple for that reason.
#
# histall is still the right option: it enables f_sig1, f_sig2, f_sigP and
# f_trsig, read inside the write_history gate, AND f_strintx, f_strinty and
# f_taubx, which ice_history accumulates every timestep -- so both halves of
# the scatter decision, deferred and not-deferred, are exercised.  run2day
# dumps a restart at day 2, which is what gives comparebfb something to read.
smoke          gx3     4x4        histall,run2day
smoke          gx3     4x4        histall,run2day,evp1d     smoke_gx3_4x4_histall_run2day
restart        gx1    40x4        droundrobin,medium
restart        tx1    40x4        dsectrobin,medium
restart        tx1    40x4        dsectrobin,medium,jra55do
restart        gx3     4x4        medium
restart        gx3     4x4        gx3nc,short
restart        gx3    10x4        maskhalo,medium
restart        gx3     6x2        alt01
restart        gx3     8x2        alt02
restart        gx3     4x2        alt03
restart        gx3    12x2        alt03,maskhalo,droundrobin
restart        gx3     4x4        alt04
restart        gx3     4x4        alt05,medium
restart        gx3     8x2        alt06
restart        gx3     8x2        pondsealvl
restart        gx3    16x2        snicar
restart        gx3    18x2        debug,maskhalo
restart        gx3     6x2        alt01,debug,short
restart        gx3     8x2        alt02,debug,short
restart        gx3     4x2        alt03,debug,short
smoke          gx3    12x2        alt03,debug,short,maskhalo,droundrobin
smoke          gx3     4x4        alt04,debug,short
smoke          gx3     4x4        alt05,debug,short
smoke          gx3     8x2        alt06,debug,short
smoke          gx3     8x3        alt07,debug,short
smoke          gx3     8x2        congel,debug,short
smoke          gx3    16x2        snicar,debug,short
smoke          gx3    12x2        snicartest,debug,short
smoke          gx3     10x2       debug,diag1,run5day,gx3sep2
smoke          gx3     7x2x5x29x12 diag1,bigdiag,run1day,debug
restart        gbox128 4x2        short
restart        gbox128 4x2        boxnodyn,short
restart        gbox128 4x2        boxnodyn,short,debug
restart        gbox128 2x2        boxadv,short
smoke          gbox128 2x2        boxadv,short,debug
restart        gbox128 4x4        boxrestore,medium
smoke          gbox128 4x4        boxrestore,short,debug
restart        gbox80  1x1        box2001
smoke          gbox80  1x1        boxslotcyl
smoke          gbox80  8x2        boxgauss,bclinearextrap,debug
smoke          gbox80  9x2        boxgauss,bczerogradient,restore5,debug
smoke          gbox12  1x1x12x12x1  boxchan,diag1,debug
restart        gx3     8x2        modal
smoke          gx3     8x2        bgcz,diag1,run5day
smoke          gx3     8x2        jra55do
smoke          gx3     8x2        bgczm,diag1,debug
smoke          gx3    12x2        zaero,diag1,debug
#smoke          gx3     8x1        bgcskl,diag1,debug
#smoke          gx3     4x1       bgcz,thread        smoke_gx3_8x2_bgcz
#restart        gx1     4x2        bgcsklclim,medium
restart        gx1     8x1        bgczclim,medium
restart        gx3    16x1        zaero,icdefault,snwitdrdg,snwgrain
smoke          gx1    24x1        medium,run90day,yi2008
smoke          gx1    24x1        medium,run90day,yi2008,jra55do
smoke          gx3     8x1        medium,run90day,yi2008
restart        gx1    24x1        short
restart        gx1    16x2        seabedLKD,gx1apr,short,debug
restart        gx1    15x2        seabedprob
restart        gx1    32x1        gx1prod
smoke          gx3     4x2        fsd1,diag24,run5day,debug
smoke          gx3     8x2        fsd12,diag24,run5day
restart        gx3     4x2        fsd12,debug,short
smoke          gx3     8x2        fsd12ww3,diag24,run1day
smoke          gx3     4x1        isotope,debug
restart        gx3     8x2        isotope
smoke          gx3     4x1        snwitdrdg,snwgrain,icdefault,debug
smoke          gx3     4x1        snw30percent,icdefault,debug
restart        gx3     8x2        snwitdrdg,icdefault,snwgrain
restart        gx3     4x4        gx3ncarbulk,iobinary,medium
restart        gx3     4x4        cdf64,histall,precision8,medium
smoke          gx3    30x1        bgcz,histall
smoke          gx3    14x2        fsd12,histall
smoke          gx3     4x1        dynpicard
# VP coverage.  dynpicard had one smoke test here and one decomp test, and
# dynanderson -- a whole branch of the nonlinear solver, plus
# use_mean_vrel = .false. -- had none anywhere.  The restart cases are
# self-comparing (restart vs continuous), so they check state handling even
# without a baseline.  The gbox80 case is cheap and deterministic.
restart        gx3     4x2        dynpicard,diag1
smoke          gx3     4x1        dynanderson
restart        gx3     4x2        dynanderson,diag1
smoke          gbox80  4x2        boxopen,kmtislands,boxforcee,run1day,dynpicard
restart        gx3     8x2        gx3ncarbulk,debug
restart        gx3     4x4        diag1,gx3ncarbulk,short
smoke          gx3     4x1        calcdragio
restart        gx3     4x2        atmbndyconstant
restart        gx3     4x2        atmbndymixed
smoke          gx3    12x2        diag1,run5day,restaicetest,debug
