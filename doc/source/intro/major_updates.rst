:tocdepth: 3

.. _updates:


Major CICE updates
============================================

This model release is CICE version 6.0.0alpha.

~~~~~~~~~
CICE V5.1
~~~~~~~~~

- include ice velocity in atm-ice coupling updates (e.g. stress) for high-frequency coupling
- allow a variable coefficient for the ice-ocean heat flux
- several new namelist options improve flexibility, especially for coupled model configurations:
   - ice-ocean heat flux
   - 'virtual' or 'real' topo melt pond water
   - ocean freezing temperature
   - high-frequency coupling
   - coupling and computational grids may be different
   - and more
- additional software enhancements improve flexibility and compatibility with CESM, Hadley Centre, and U.S. Navy coupled models
- new diagnostics and updated documentation
- various bug fixes 

~~~~~~~~~
CICE V5.0
~~~~~~~~~

- A method for prognosing sea ice salinity, including improved snow-ice formation
- Two new explicit melt pond parameterizations (topo and level-ice)
- Sea ice biogeochemistry
- Elastic-Anisotropic-Plastic rheology
- Improved parameterization of form drag 
- The "revised EVP" under-damping approach
- Gracefully handles the case when an internal layer melts completely
- Gregorian calendar with leap years
- Reduced memory and floating-point operations for tracer calculations
- Ice and snow enthalpy defined as tracers
- New history variables for melt ponds, ridging diagnostics, biogeochemistry and more
- Read/write variables on the extended grid, including ghost (halo) cells
- Parallelism option via OpenMP threads
- Improved parallel efficiency through halo masks and new block distributions
- Parallel I/O option via the PIO library
- Restarts in binary or netCDF formats
- CPP options for categories, layers and tracers 
- Corrected bugs, particularly for nonstandard configurations.