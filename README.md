
## Overview

This repository contains files needed to run versions 6 and higher of the sea ice model CICE, which is now maintained by the CICE Consortium.  Versions prior to v6 are found in the [CICE-svn-trunk repository](https://github.com/CICE-Consortium/CICE-svn-trunk).

CICE consists of a top level driver and dynamical core plus the Icepack column physics code, which is included in CICE as a git submodule.  Because Icepack is a submodule of CICE, Icepack and CICE development are handled independently with respect to the github repositories even though development and testing may be done together. 

## Obtaining CICE

A list of the official CICE releases along with release notes is located here:
https://github.com/CICE-Consortium/CICE/releases

If you use the CICE code, please cite the version you are using with the CICE Digital Object Identifier (DOI):
DOI:10.5281/zenodo.1205674 (https://zenodo.org/record/1205674)
This DOI can be used to cite all CICE versions and the URL will default to the most recent version. However, each released version of CICE will also receive its own, unique DOI that can be used for citations as well.

CICE may be obtained in several different ways:  [not yet tested]    
1.  clone the full repository
See [Git Workflow Guidance](https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance)
2.  check out only a particular branch, version or tag
In the workflow for step 1 above, substitute
git clone -b branch_name https://github.com/CICE-Consortium/CICE.git local_directory_name
or use svn
svn co https://github.com/CICE-Consortium/CICE/branch_name
where "branch name" can also be a version name
3.  download a tarball for a particular version from the git releases: https://github.com/CICE-Consortium/CICE/releases

If you expect to make any changes to the code, we recommend that you first fork both the CICE and Icepack repositories.  Basic instructions for working with CICE and Icepack are found in the [Git Workflow Guidance](https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance), linked from the wikis in the primary code repositories
https://github.com/CICE-Consortium/CICE/wiki
https://github.com/CICE-Consortium/Icepack/wiki

## Documentation

Detailed and searchable online documentation of CICE can be found at https://readthedocs.org/projects/cice-consortium-cice/ .

This site has the most up-to-date [HTML](http://cice-consortium-cice.readthedocs.io/en/master/) and [PDF](https://media.readthedocs.org/pdf/cice-consortium-cice/master/cice-consortium-cice.pdf) living documentation from the master branch of the CICE-Consortium repository that will be updated regularly with code development.

This site also has static documentation from each CICE release.

More information about CICE documentation can be found on the [CICE Documentation Wiki page](https://github.com/CICE-Consortium/CICE/wiki/CICE-Documentation).

## More Information

The [CICE wiki](https://github.com/CICE-Consortium/CICE/wiki) page contains links to additional information, e.g.
- larger files such as the gx1 grid, land mask, and forcing files
- testing data

The [Test-Results wiki](https://github.com/CICE-Consortium/Test-Results/wiki) has test results for both CICE and Icepack.

The [About-Us repository](https://github.com/CICE-Consortium/About-Us) includes background and supporting information about the CICE Consortium, including how to interact with it.   

See also our [FAQ](https://github.com/CICE-Consortium/About-Us/wiki/FAQ).



