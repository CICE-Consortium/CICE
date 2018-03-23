:tocdepth: 3 

.. _about_dev:

About Development
==================

The Icepack model consists of three different parts, the column physics
code, the icepack driver, and the scripts.  Development of each of these
pieces will be described below separately.

Subroutine calls and other linkages into Icepack from the host model should only
need to access the **icepack\_intfc\*.F90** interface modules within the 
``columnphysics/`` directory.  
The Icepack driver in the ``configuration/driver/`` directory is based on the CICE
model and provides an example of the sea ice host model capabilities needed for inclusion
of Icepack.  In particular, host models will need to include code equivalent to that
in the modules **icedrv\_\*_column.F90**.  Calls into the Icepack interface routines
are primarily from **icedrv\_step\_mod.F90** but there are others (search the driver code
for ``intfc``).

Guiding principles for the creation of Icepack include the following: 
  - The column physics modules shall be independent of all sea ice model infrastructural
    elements that may vary from model to model.  Examples include input/output, timers,
    references to CPUs or computational tasks, initialization other than that necessary for
    strictly physical reasons, and anything related to a horizontal grid.
  - The column physics modules shall not call or reference any routines or code that 
    reside outside of the **columnphysics/** directory.
  - Any capabilities required by a host sea ice model (e.g. calendar variables, tracer 
    flags, diagnostics) shall be implemented in the driver and passed into or out of the 
    column physics modules via array arguments.
