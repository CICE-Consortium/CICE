:tocdepth: 3

.. _dev_dynamics:


Dynamics
============================

The CICE **cicecore/** directory consists of the non icepack source code.  Within that
directory there are the following subdirectories

**cicecore/cicedyn/analysis** contains higher level history and diagnostic routines.

**cicecore/cicedyn/dynamics** contains all the dynamical evp, eap, and transport routines.

**cicecore/cicedyn/general** contains routines associated with forcing, flux calculation,
initialization, and model timestepping.

**cicecore/cicedyn/infrastructure** contains most of the low-level infrastructure associated
with communication (halo updates, gather, scatter, global sums, etc) and I/O reading and writing
binary and netcdf files.

**cicecore/drivers/** contains subdirectories that support stand-alone drivers and other high level
coupling layers.

**cicecore/shared/** contains some basic methods related to grid decomposition, time managers, constants, kinds, and restart capabilities.


Dynamical Solvers
--------------------

The dynamics solvers are found in **cicecore/cicedyn/dynamics/**.  A couple of different solvers are
available including EVP, EAP and VP.  The dynamics solver is specified in namelist with the
``kdyn`` variable.  ``kdyn=1`` is evp, ``kdyn=2`` is eap, ``kdyn=3`` is VP.

Two alternative implementations of EVP are included. The first alternative is the Revised EVP, triggered when the ``revised_evp`` is set to true. The second alternative is the 1d EVP solver triggered when the ``evp_algorithm`` is set to ``shared_mem_1d`` as oppose to the default setting of ``evp_standard_2d``. The solutions with ``evp_algorithm`` set to ``standard_2d`` or ``shared_mem_1d`` will
not be bit-for-bit identical when compared to each other. The reason for this is floating point round off errors that occur unless strict compiler flags are used. ``evp_algorithm=shared_mem_1d`` is primarily built for OpenMP. If MPI domain splitting is used then the solver will only run on the master processor. ``evp_algorithm=shared_mem_1d`` is not supported
with the tripole grid.


Transport
-----------------

The transport (advection) methods are found in **cicecore/cicedyn/dynamics/**.  Two methods are supported,
upwind and remap.  These are set in namelist via the ``advection`` variable.
Transport can be disabled with the ``ktransport`` namelist variable.


