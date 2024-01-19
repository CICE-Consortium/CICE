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


Coding Standard
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Overall, CICE code should be implemented as follows,

  * Adhere to the current coding and naming conventions

  * Write readable code.  Use meaningful variable names; indent 2 or 3 spaces for loops and conditionals; vertically align similar elements where it makes sense, and provide concise comments throughout the code.

  * Declare common parameters in a shared module.  Do not hardwire the same parameter in the code in multiple places.

  * Maintain bit-for-bit output for the default configuration (to the extent possible).  Use namelist options to add new features.

  * Maintain global conservation of heat, water, salt

  * Use of C preprocessor (CPP) directives should be minimized and only used for build dependent modifications such as use of netcdf (or other "optional" libraries) or for various Fortran features that may not be supported by some compilers. Use namelist to support run-time code options. CPPs should be all caps.

  * All modules should have the following set at the top

    .. code-block:: fortran

       implicit none
       private

    Any public module interfaces or data should be explicitly specified

  * All subroutines and functions should define the ``subname`` character parameter statement to match the interface name like

    .. code-block:: fortran

       character(len=*),parameter :: subname='(advance_timestep)'

  * Public Icepack interfaces should be accessed thru the ``icepack_intfc`` module like

    .. code-block:: fortran

       use icepack_intfc, only: icepack_init_parameters

  * Icepack does not write to output or abort, it provides methods to access those features.  After each call to Icepack, **icepack_warnings_flush** should be called to flush Icepack output to the CICE log file and **icepack_warnings_aborted** should be check to abort on an Icepack error as follows,

    .. code-block:: fortran

       call icepack_physics()
       call icepack_warnings_flush(nu_diag)
       if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)

  * Use ``ice_check_nc`` or ``ice_pio_check`` after netcdf or pio calls to check for return errors.

  * Use subroutine ``abort_ice`` to abort the model run. Do not use stop or MPI_ABORT.  Use optional arguments (file=__FILE__, line=__LINE__) in calls to ``abort_ice`` to improve debugging

  * Write output to stdout from the master task only unless the output is associated with an abort call.  Write to unit ``nu_diag`` following the current standard.  Do not use units 5 or 6.  Do not use the print statement.

  * Use of new Fortran features or external libraries need to be balanced against usability and the desire to compile on as many machines and compilers as possible.  Developers are encouraged to contact the Consortium as early as possible to discuss requirements and implementation in this case.

