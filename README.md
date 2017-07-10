*This repository is under construction*

## Overview

This repository contains files needed to run versions 6 and higher of the sea ice model CICE, which is now maintained by the CICE Consortium.  Versions prior to v6 are found in the CICE-svn-trunk repository    
https://github.com/CICE-Consortium/CICE-svn-trunk 

CICE consists of a top level driver and dynamical core plus the Icepack column physics code, which is included in CICE as a git submodule.  Because Icepack is a submodule of CICE, Icepack and CICE development are handled independently with respect to the github repositories even though development and testing may be done together. 

## Obtaining CICE

If you expect to make any changes to the code, we recommend that you first fork both the CICE and Icepack repositories.  Basic instructions for working with CICE and Icepack are found in the Git and Workflow Guide, linked from the wikis in the primary repositories    
https://github.com/CICE-Consortium/CICE    
https://github.com/CICE-Consortium/Icepack

CICE may be obtained in several different ways:  [not yet tested]    
1.  clone the full repository    
See Git and Workflow Guide    
2.  check out only a particular branch, version or tag    
In the workflow for step 1 above, substitute    
  [check this] git clone -b branch_name --single-branch --recursive https://github.com/CICE-Consortium/CICE.git local_directory_name
or use svn    
  svn co https://github.com/CICE-Consortium/CICE/branch_name    
where “branch name” can also be a version name    
3.  download a tarball for a particular version    
[how]

## More information

"Quick Start" instructions are available in README_v6, and instructions for setting up standard tests (e.g. regression, restart) are in README_test.  

The doc directory contains scientific documentation.

 [check this]   The wiki pages for each repository contain links to additional information, e.g.    
- complete documentation 
- larger files such as the gx1 grid, land mask, and forcing files
- testing data
- test results 

The "About-Us" repository includes background and supporting information about the CICE Consortium, including how to interact with it.    
https://github.com/CICE-Consortium/About-Us
