:tocdepth: 3


.. _quickstart:

Quick Start
===========

Clone the model from the CICE-Consortium repository::

  git clone --recurse-submodules https://github.com/CICE-Consortium/CICE

Instructions for working with Git and GitHub with CICE (and Icepack) can be
found in the `CICE Git Workflow Guide <https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guide>`_.

You will probably have to download some input data, see the `CICE wiki <https://github.com/cice-consortium/CICE/wiki>`_ or :ref:`force`.

Software requirements are noted in this :ref:`software` section.

Porting information can be found in the :ref:`porting` section.  A special porting section for personal computers 
is in the :ref:`laptops` section.

From your main CICE directory, execute::

  ./cice.setup -c ~/mycase1 -g gx3 -m MACHINE -e ENVIRONMENT -s diag1,thread -p 8x1
  cd ~/mycase1
  ./cice.build
  ./cice.submit

where ``MACHINE`` is the machine name and ``ENVIRONMENT`` is the compiler environment,
for example ``-m derecho -e intel``.
The supported combinations are the **env.[machine]_[environment]** files in
**configuration/scripts/machines/**, and ``./cice.setup --help`` lists the setup options.
See :ref:`scripts` for more information about how to use the cice.setup and cice.submit scripts.

To run CICE on a laptop or personal computer, use the conda configuration
(``-m conda -e macos`` or ``-m conda -e linux``), see :ref:`laptops`.
To run on another HPC system, the model has to be ported first, see :ref:`porting`.
In both cases the input data must also be downloaded, see :ref:`force`.

Please cite any use of the CICE code. More information can be found at :ref:`citing`.

