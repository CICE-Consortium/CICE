:tocdepth: 3 

.. _dev_icepack:

Icepack
==================

The CICE model calls the Icepack columnphysics source code.  The Icepack model is documented
separately, see https://github.com/CICE-Consortium/Icepack.

More specifically, the CICE model uses methods defined in **icepack_intfc.F90**.  It uses 
the init, query, and write methods to set, get, and document Icepack values.  And it follows
the icepack_warnings methodology where icepack_warnings_aborted is checked and
icepack_warnings_print is called after every call to an Icepack method.  It does not directly
"use" Icepack data and access Icepack data only thru interfaces.



