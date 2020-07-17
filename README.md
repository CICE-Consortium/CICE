[![Build Status](https://travis-ci.org/CICE-Consortium/CICE.svg?branch=master)](https://travis-ci.org/CICE-Consortium/CICE)
[![Documentation Status](https://readthedocs.org/projects/cice-consortium-cice/badge/?version=master)](http://cice-consortium-cice.readthedocs.io/en/master/?badge=master)
[![lcov](https://img.shields.io/endpoint?url=https://apcraig.github.io/coverage.json)](https://apcraig.github.io)

<!--- [![codecov](https://codecov.io/gh/apcraig/Test_CICE_Icepack/branch/master/graph/badge.svg)](https://codecov.io/gh/apcraig/Test_CICE_Icepack) --->

## The CICE Consortium sea-ice model
CICE is a computationally efficient model for simulating the growth, melting, and movement of polar sea ice. Designed as one component of coupled atmosphere-ocean-land-ice global climate models, today’s CICE model is the outcome of more than two decades of community collaboration in building a sea ice model suitable for multiple uses including process studies, operational forecasting, and climate simulation.


This repository contains the files and code needed to run the CICE sea ice numerical model starting with version 6. CICE is maintained by the CICE Consortium. 
Versions prior to v6 are found in the [CICE-svn-trunk repository](https://github.com/CICE-Consortium/CICE-svn-trunk).

CICE consists of a top level driver and dynamical core plus the [Icepack][icepack] column physics code], which is included in CICE as a Git submodule.  Because Icepack is a submodule of CICE, Icepack and CICE development are handled independently with respect to the GitHub repositories even though development and testing may be done together.  

[icepack]: https://github.com/CICE-Consortium/Icepack

The first point of contact with the CICE Consortium is the Consortium Community [Forum][forum]. 
This forum is monitored by Consortium members and also opened to the whole community.
Please do not use our issue tracker for general support questions.

[forum]: https://xenforo.cgd.ucar.edu/cesm/forums/cice-consortium.146/

If you expect to make any changes to the code, we recommend that you first fork both the CICE and Icepack repositories. 
In order to incorporate your developments into the Consortium code it is imperative you follow the guidance for Pull Requests and requisite testing.
Head over to our [Contributing][contributing] guide to learn more about how you can help improve CICE.

[contributing]: https://github.com/CICE-Consortium/About-Us/wiki/Contributing

## Useful links
* **CICE wiki**: https://github.com/CICE-Consortium/CICE/wiki

   Information about the CICE model

* **CICE Release Table**: https://github.com/CICE-Consortium/CICE/wiki/CICE-Release-Table

   Numbered CICE releases since version 6 with associated documentation and DOIs. 
   
* **Consortium Community Forum**: https://xenforo.cgd.ucar.edu/cesm/forums/cice-consortium.146/

   First point of contact for discussing model development including bugs, diagnostics, and future directions.   

* **Resource Index**: https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index

   List of resources for information about the Consortium and its repositories as well as model documentation, testing, and development.

## License
See our [License](LICENSE.pdf) and [Distribution Policy](DistributionPolicy.pdf).
