:tocdepth: 3


.. _quickstart:

Quick Start
===========

Download the model from the CICE-Consortium repository, 
    https://github.com/CICE-Consortium/CICE

Instructions for working in github with CICE (and Icepack) can be
found in the `CICE Git and Workflow Guide <https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance>`_.

You will probably have to download some inputdata, see the `CICE wiki <https://github.com/cice-consortium/CICE/wiki>`_ or :ref:`force`.

From your main CICE directory, execute::

  ./cice.setup -c ~/mycase1 -g gx3 -m testmachine -s diag1,thread -p 8x1
  cd ~/mycase1
  ./cice.build
  ./cice.submit


``testmachine`` is a generic machine name included with the cice scripts.
The local machine name will have to be substituted for ``testmachine`` and
there are working ports for several different machines.  However, it may be necessary
to port the model to a new machine.  See :ref:`porting` for 
more information about how to port. See :ref:`scripts` for more information about 
how to use the cice.setup and cice.submit scripts.

Please cite any use of the CICE code. More information can be found at :ref:`citing`.

