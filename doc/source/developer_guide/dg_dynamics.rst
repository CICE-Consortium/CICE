:tocdepth: 3

.. _dev_dynamics:


Dynamics and Infrastructure Implementation
================================================

The CICE **cicecore/** directory consists of the non icepack source code.  Within that 
directory there are the following subdirectories

**cicedynB/analysis** contains higher level history and diagnostic routines.

**cicedynB/dynamics** contains all the dynamical evp, eap, and transport routines.

**cicedynB/general** contains routines associated with forcing, flux calculation,
initialization, and model timestepping.

**cicedynB/infrastructure** contains most of the low-level infrastructure associated
with communication (halo updates, gather, scatter, global sums, etc) and I/O reading and writing
binary and netcdf files.

**drivers/** contains subdirectories that support stand-alone drivers and other high level
coupling layers.

**shared/** contains some basic methods related to grid decomposition, time managers, constants,
kinds, and restart capabilities.


Dynamics
~~~~~~~~~~~~~~

Dyanamical Solvers
************************

EVP and EAP plus revised EVP


Transport
**************

remap and other stuff?


Infrastructure
~~~~~~~~~~~~~~~~~~~~

Kinds
*********

Constants
*************

Static Array Allocation
**************************

Time Manager
****************

Communication
********************

I/O
***********

