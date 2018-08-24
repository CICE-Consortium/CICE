[![Build Status](https://travis-ci.org/CICE-Consortium/CICE.svg?branch=master)](https://travis-ci.org/CICE-Consortium/CICE)
[![Documentation Status](https://readthedocs.org/projects/cice-consortium-cice/badge/?version=master)](http://cice-consortium-cice.readthedocs.io/en/master/?badge=master)

## Overview
This repository contains the files and code needed to run the CICE sea ice numerical model starting with version 6. CICE is maintained by the CICE Consortium. Versions prior to v6 are found in the [CICE-svn-trunk repository](https://github.com/CICE-Consortium/CICE-svn-trunk).

CICE consists of a top level driver and dynamical core plus the Icepack column physics code, which is included in CICE as a git submodule.  Because Icepack is a submodule of CICE, Icepack and CICE development are handled independently with respect to the github repositories even though development and testing may be done together. 

If you expect to make any changes to the code, we recommend that you first fork both the CICE and Icepack repositories. Basic instructions for working with CICE and Icepack are found in the Git Workflow Guidance, linked from the Resource Index (below). In order to incorporate your developments into the Consortium code it is
imperative you follow the guidance for Pull Requests and requisite testing.

## Useful links
* **CICE wiki**: https://github.com/CICE-Consortium/CICE/wiki

   Information about the CICE model

* **CICE Version Index**: https://github.com/CICE-Consortium/CICE/wiki/CICE-Version-Index

   Numbered CICE releases since version 6 with associated documentation and DOIs. 

* **Resource Index**: https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index

   List of resources for information about the Consortium and its repositories as well as model documentation, testing, and development.
