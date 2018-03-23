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

This model release is CICE version 6.0.0alpha, available from https://github.com/CICE-Consortium/CICE.
It updates CICE5.1, and there are considerable changes to the code.  

This document uses the following text conventions:
Variable names used in the code are ``typewritten``.
Subroutine names are given in *italic*.
File and directory names are in **boldface**.
A comprehensive :ref:`index`, including glossary of symbols with many of their values, appears
at the end of this guide.