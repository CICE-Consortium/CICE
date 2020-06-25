:tocdepth: 3

.. _about:

About CICE
=============
CICE is a computationally efficient model for simulating the growth, 
melting, and movement of polar sea ice. Designed as one component of 
coupled atmosphere-ocean-land-ice global climate models, today’s CICE 
model is the outcome of more than two decades of effort led by 
scientists at Los Alamos National Laboratory. The current version of 
the model has been enhanced greatly through collaborations with members 
of the community.

CICE has several interacting components: a model of ice dynamics, which 
predicts the velocity field of the ice pack based on a model of the 
material strength of the ice; a transport model that describes advection 
of the areal concentration, ice volumes and other state variables; and a 
vertical physics package, called “Icepack”, which includes mechanical, 
thermodynamic, and biogeochemical models to compute thickness changes 
and the internal evolution of the hydrological ice-brine ecosystem. When 
coupled with other earth system model components, routines external to the 
CICE model prepare and execute data exchanges with an external “flux coupler”.

Icepack is implemented in CICE as a git submodule, and it is documented at 
https://cice-consortium-icepack.readthedocs.io/en/master/index.html. 
Development and testing of CICE and Icepack may be done together,
but the repositories are independent.
This document describes the remainder of the CICE model. The CICE code is 
available from https://github.com/CICE-Consortium/CICE.

The standard standalone CICE test configuration uses a 3 degree grid with 
atmospheric data from 1997, available at
https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data.
A 1-degree configuration and data are also available, along with some idealized 
configurations. The data files are designed only for testing the code, not 
for use in production runs or as observational data. Please do not publish 
results based on these data sets.

The CICE model can run serially or in parallel, and the CICE software package 
includes tests for various configurations. MPI is used for message passing 
between processors, and OpenMP threading is available.

Major changes with each CICE release (https://github.com/CICE-Consortium/CICE/releases) 
will be detailed with the included release notes. Enhancements and bug fixes made to 
CICE since the last numbered release can be found on the CICE wiki
(https://github.com/CICE-Consortium/CICE/wiki/CICE-Recent-changes).
**Please cite any use of the CICE code.** More information can be found at :ref:`citing`. 

This document uses the following text conventions: Variable names used in 
the code are ``typewritten``. Subroutine names are given in *italic*. File 
and directory names are in **boldface**. A comprehensive :ref:`index`, 
including glossary of symbols with many of their values, appears at the 
end of this guide.
