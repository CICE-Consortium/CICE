:tocdepth: 3

.. _about:

About CICE
=============

The Los Alamos sea ice model (CICE) is the result of an effort to
develop a computationally efficient sea ice component for a fully
coupled atmosphere--land global climate model. It was
designed to be compatible with the Parallel Ocean Program
(POP), an ocean circulation model developed at 
Los Alamos National Laboratory for use on massively parallel computers
:cite:`Smith92,Dukowicz93,Dukowicz94`. The current version of the
model has been enhanced greatly through collaborations with members of
the community.

CICE has several interacting components: a thermodynamic model that
computes local growth rates of snow and ice due to vertical conductive, 
radiative and turbulent fluxes, along with snowfall and melt ponds; a 
model of ice dynamics, which predicts the velocity field of the ice pack 
based on a model of the material strength of the ice; a transport model 
that describes advection of the areal concentration, ice volumes and other 
state variables; and a ridging parameterization that transfers ice among
thickness categories based on energetic balances and rates of strain. 
When coupled with other earth system model components, routines external
to the CICE model prepare and execute data exchanges with an external
"flux coupler".

Details about this model release and a list of major changes are found 
in :ref:`updates` and the model code
is available from https://github.com/CICE-Consortium/CICE. 

Please cite any use of the CICE code. More information can be found at :ref:`citing`.

The standard standalone CICE test configuration uses a 3 degree grid with 
atmospheric data from 1997, available as detailed on the 
`GitHub wiki <https://github.com/CICE-Consortium/CICE/wiki>`_. 
A 1 degree configuration
and data are also available, along with some idealized configurations. The
data files are designed only for testing the code, not for use in production
runs or as observational data. **Please do not publish results based on these
data sets.** 

The CICE model can run serially or in parallel, and the CICE software package
includes tests for various configurations. MPI is used for message passing
between processors, and OpenMP threading is available.

This document uses the following text conventions:
Variable names used in the code are ``typewritten``.
Subroutine names are given in *italic*.
File and directory names are in **boldface**.
A comprehensive :ref:`index`, including glossary of symbols with many of their values, appears
at the end of this guide.
