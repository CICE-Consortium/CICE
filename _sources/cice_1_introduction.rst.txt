Introduction - CICE5
============================================

The Los Alamos sea ice model (CICE) is the result of an effort to
develop a computationally efficient sea ice component for a fully
coupled atmosphere--land global climate model. It was
designed to be compatible with the Parallel Ocean Program
(POP), an ocean circulation model developed at 
Los Alamos National Laboratory for use on massively parallel computers
:cite:`SDM92,DSM93,DSM94`. The current version of the
model has been enhanced greatly through collaborations with members of
the community.

CICE has several interacting components: a thermodynamic model that
computes local growth rates of snow and ice due to vertical conductive, 
radiative and turbulent fluxes, along with snowfall; a model of ice 
dynamics, which predicts the velocity field of the ice pack based on 
a model of the material strength of the ice; a transport model that 
describes advection of the areal concentration, ice volumes and other 
state variables; and a ridging parameterization that transfers ice among
thickness categories based on energetic balances and 
rates of strain.External routines would prepare and execute data exchanges with an
external "flux coupler," which then passes the data to other climate
model components such as POP.

This model release is CICE version 5.1, available from http://oceans11.lanl.gov/trac/CICE/wiki.
It updates CICE5.0, which was released in September 2013. With so many new parameterizations,
we must admit that all combinations have not been tested.  Also, different parameterizations for 
various sub-processes (e.g., snow infiltration by melt or sea water) have been introduced as part 
of the new code options, which need to be unified or pulled out as separate options.  

This document uses the following text conventions:
Variable names used in the code are ``typewritten``.
Subroutine names are given in *italic*.
File and directory names are in **boldface**.
A comprehensive :ref:`index`, including glossary of symbols with many of their values, appears
at the end of this guide.

======================
Major updates in V5.1
======================
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

======================
Major updates in V5.0
======================
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

