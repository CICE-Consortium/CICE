:tocdepth: 3 

.. _dev_about:

About Development
==================

The CICE model consists of four different parts, the CICE dynamics and supporting infrastructure, 
the CICE driver code, the Icepack column physics code, and the scripts.  Development of each of these
pieces is described separately.

Guiding principles for the creation of CICE include the following: 
  - CICE can be run in stand-alone or coupled modes.  A top layer driver, coupling layer,
    or model cap can be used to drive the CICE model.
  - The Icepack column physics modules are independent, consist of methods that operate
    on individual gridcells, and contain no underlying infrastructure.  CICE must call
    into Icepack using interfaces and approaches specified by Icepack.


Git workflow and Pull Requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is extensive Information for Developers documentation available.  See https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index#information-for-developers for information on:
  - Contributing to model development
  - Software development practices guide
  - git Workflow Guide - including extensive information about the Pull Request process and requirements
  - Documentation Workflow Guide
