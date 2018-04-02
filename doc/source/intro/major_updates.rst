:tocdepth: 3

.. _updates:


Major CICE updates
============================================

This model release is CICE version 6.0.0.alpha.

Please cite any use of the CICE code. More information can be found at :ref:`citing`.

~~~~~~~~~~~~~~~~~
CICE V6.0.0.alpha
~~~~~~~~~~~~~~~~~
Major changes:

- A new fast-ice parameterization
- Full vertical biogeochemistry
- Independent column physics package Icepack implemented as a git submodule
- A flexible, extensible, robust interface between the column physics modules and the driver
- A warning package that captures diagnostic and error information from within the column physics, for printing by the driver
- Restructured code and forcing data directories
- An entirely new scripting system
- A comprehensive test suite of various configuration options, with quality control and compliance tests
- Automated testing using Travis CI
- Automated test reporting organized by hash, version, machine and branch, for both the primary Consortium repository and user forks
- Online documentation
- See also updates in Icepack releases and recent changes

Enhancements:

- Change use of ULAT to TLAT to determine what latitudes initial ice is present in set_state_var [r970]
- Add 4d fields to history (categories, vertical ice) r1076
- Update PIO; Universal large file support [r1094]
- Remove calendar_type from namelist options and initialize it based on the namelist flag use_leap_years. [r1098]
- Add fbot to history output [r1107]
- Add shortwave diagnostics [r1108]
- Modifications to enable ocean and ice biogeochemical coupling [r1111, r1200]
- Remove the computational overhead of coupling BGC when it is not being used [r1123]
- Change reference to char_len in stop_label [r1143]
- Add grounding scheme and tensile strength #52
- Add new namelist options for dynamics parameters #52
- Update Icepack version in CICE (Icepack v1.0.0 #81)
- Modifications to stress diagnostics, including principal stress normalization and internal pressure #99

Bug fixes:

- Properly read and rotate ocean currents from 3D gx1 netcdf data r959
- Correct diagnostic output 'avg salinity' [r1022]
- Bug fix for padded domains. r1031
- Use VGRD instead of VGRDi for 3D [r1037]
- change shortwave calculation to depend on the net shortwave sum instead of cosine of the zenith angle (not BFB: in addition to the different shortwave calculation, albedo output in history is different). r1076
- Correct available history fields. [r1082]
- Fix coupled restart bug; initialize coszen; adjust calendar_type implementation [r1094]
- Port miscellaneous changes from the various column package branches back to the trunk. BFB in the standard configuration, but the initializations and conditional changes for coszen could change the answers in other configurations. Also the flux calculation change in ice_therm_itd.F90 could change the answers in coupled simulations. 1102
- Ensure fractions of snow, ponds and bare ice add to one r1120
- Zero out thin-pond fraction for radiation in cesm, topo pond schemes (not BFB), and set albedo=1 where/when there is no incoming shortwave (changes the average-albedo diagnostic), and fix thin (cesm) ponds overlapping snow. [r1126, r1132]
- Fix padding when using the extended-grid functionality, to prevent arrays out of bounds. [r1128]
- Change dynamics halo update mask from icetmask to iceumask (fixes occasional exact restart problem and error in halo update) [r1133]
- Add surface flooding and surface runoff terms which increase with open water area in surface condition for update_hbrine, z_salinity, z_biogeochemistry [r1161]
- Set all tracer values to c0 over land after initialization #16
- Remove OpenMP directives for loops that do not appear to be thread safe #25
- Remove iblk from timer starts #98
