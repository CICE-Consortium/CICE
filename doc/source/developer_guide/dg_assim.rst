:tocdepth: 3

.. _dataassimilation:

Data Assimilation
======================

Data assimilation (DA) is the scientific process of combining external
data with numerical model forecasts.  There are several ways this can
be done including by adjusting the model initial conditions (internally
or externally) or adjusting the model solution as it evolves in time.
Various data assimilation options are being introduced in CICE and are
described below.

.. _restartmod:

Data Assimilation on restart
------------------------------------

The namelist variable, ``restart_mod``, specifies the restart DA mode.  
By default, this namelist value is set to ``none`` which disables the feature.  
The current active options are ``adjust_aice`` and ``adjust_aice_test``.

With ``adjust_aice`` and ``adjust_aice_test``, the category averaged aice 
value is modified at restart to specified values using the method implemented in 
**cicecore/cicedyn/infrastructure/ice_restart_driver.F90** subroutine
**direct_adjust_aice**.  This method adjusts aice, vice, vsno, qice, and
sice in all categories to be consistent with the category average aice
specified.  It also adjusts several thermodynamic variables such as 
temperature and salinity (see :cite:`Posey15`).  
``adjust_aice`` reads in a sea ice concentration
field from an external file.  The field is currently hardwired to 'sic' and the
file is currently hardwired to 'sic.nc'.  The field must be on the model grid.
``adjust_aice_test`` modifies the
aice field read on restart internally.  The current implementation rounds
the aice values read at restart to the nearest 1/100th.  This mode exists
primarily to test the feature.
