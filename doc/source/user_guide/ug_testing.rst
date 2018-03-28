:tocdepth: 3

Testing CICE
============

Version 6, August 2017
This documents how to use the testing features developed for the 
CICE Consortium CICE sea ice model.

.. _basic:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Individual tests and test suites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CICE scripts support both setup of individual tests as well as test suites.  Individual
tests are run from the command line like

  > cice.setup -t smoke -m wolf -g gx3 -p 8x2 -s diag1,run5day -testid myid

where -m designates a specific machine.  Test suites are multiple tests that are specified in 
an input file and are started on the command line like

  > cice.setup -ts base_suite -m wolf -testid myid

cice.setup with -t or -ts require a testid to uniquely name test directories.  The format
of the case directory name for a test will always be 
${machine}_${test}_${grid}_${pes}_${soptions}.${testid}

To build and run a test, the process is the same as a case,
  cd into the test directory,
  
  run cice.build
  
  run cice.submit

The test results will be generated in a local file called "test_output".

When running a test suite, the cice.setup command line automatically generates all the tests
under a directory names ${test_suite}.${testid}.  It then automatically builds and submits all
tests.  When the tests are complete, run the results.csh script to see the results from all the
tests.

Tests are defined under configuration/scripts/tests.  The tests currently supported are:
  smoke   - Runs the model for default length.  The length and options can
            be set with the -s commmand line option.  The test passes if the
            model completes successfully.
  restart - Runs the model for 10 days, writing a restart file at day 5 and
            again at day 10.  Runs the model a second time starting from the
            day 5 restart and writing a restart at day 10 of the model run.
            The test passes if both the 10 day and 5 day restart run complete and
            if the restart files at day 10 from both runs are bit-for-bit identical.

Please run './cice.setup -h' for additional details.

.. _additional:

~~~~~~~~~~~~~~~~~~~~~~~~~~
Additional testing options
~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several additional options on the cice.setup command line for testing that
provide the ability to regression test and compare tests to each other.

  -bd defines a baseline directory where tests can be stored for regression testing
  
  -bg defines a version name that where the current tests can be saved for regression testing
  
  -bc defines a version name that the current tests should be compared to for regression testing
  
  -td provides a way to compare tests with each other

To use -bg,
  > cice.setup -ts base_suite -m wolf -testid v1 -bg version1 -bd $SCRATCH/CICE_BASELINES
    will copy all the results from the test suite to $SCRATCH/CICE_BASELINES/version1.

To use -bc,
  > cice.setup -ts base_suite -m wolf -testid v2 -bc version1 -bd $SCRATCH/CICE_BASELINES
    will compare all the results from this test suite to results saved before in $SCRATCH/CICE_BASELINES/version1.

-bc and -bg can be combined,
  >cice.setup -ts base_suite -m wolf -testid v2 -bg version2 -bc version1 -bd $SCRATCH/CICE_BASELINES
   will save the current results to $SCRATCH/CICE_BASELINES/version2 and compare the current results to
   results save before in $SCRATCH/CICE_BASELINES/version1.

-bg, -bc, and -bd are used for regression testing.  There is a default -bd on each machine.

-td allows a user to compare one test result to another.  For instance,
  > cice.setup -t smoke -m wolf -g gx3 -p 8x2 -s run5day -testid t01
  > cice.setup -t smoke -m wolf -g gx3 -p 4x2 -s run5day -testid t01 -td smoke_gx3_8x2_run5day

An additional check will be done for the second test (because of the -td argument), and it will compare
the output from the first test "smoke_gx3_8x2_run5day" to the output from it's test "smoke_gx3_4x2_run5day"
and generate a result for that.  It's important that the first test complete before the second test is 
done.  Also, the -td option works only if the testid and the machine are the same for the baseline
run and the current run.

.. _format:

~~~~~~~~~~~~~~~~~
Test suite format
~~~~~~~~~~~~~~~~~

The format for the test suite file is relatively simple.  It is a text file with white space delimited 
columns like,

.. _tab-test:

.. csv-table:: Table 7
   :header: "#Test", "Grid", "PEs", "Sets", "BFB-compare"
   :widths: 7, 7, 7, 15, 15

   "smoke", "gx3", "8x2", "diag1,run5day", ""
   "smoke", "gx3", "8x2", "diag24,run1year,medium", ""
   "smoke", "gx3", "4x1", "debug,diag1,run5day", ""
   "smoke", "gx3", "8x2", "debug,diag1,run5day", ""
   "smoke", "gx3", "4x2", "diag1,run5day", "smoke_gx3_8x2_diag1_run5day"
   "smoke", "gx3", "4x1", "diag1,run5day,thread", "smoke_gx3_8x2_diag1_run5day"
   "smoke", "gx3", "4x1", "diag1,run5day", "smoke_gx3_4x1_diag1_run5day_thread"
   "restart", "gx3", "8x1", "", ""
   "restart", "gx3", "4x2", "debug", ""


The first column is the test name, the second the grid, the third the pe count, the fourth column is
the -s options and the fifth column is the -td argument.  The fourth and fifth columns are optional.
The argument to -ts defines which filename to choose and that argument can contain a path.  cice.setup 
will also look for the filename in configuration/scripts/tests where some preset test suites are defined.

~~~~~~~~~~~~~~~~~~~~~~~~~~
Example Tests (Quickstart)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**********************************************
To generate a baseline dataset for a test case
**********************************************

./cice.setup -t smoke -m wolf -bg cicev6.0.0 -testid t00

cd wolf_smoke_gx3_4x1.t00

./cice.build

./cice.submit

# After job finishes, check output

cat test_output

****************************************************
To run a test case and compare to a baseline dataset
****************************************************

./cice.setup -t smoke -m wolf -bc cicev6.0.0 -testid t01

cd wolf_smoke_gx3_4x1.t01

./cice.build

./cice.submit

# After job finishes, check output

cat test_output

*********************************************
To run a test suite to generate baseline data
*********************************************

./cice.setup -m wolf -ts base_suite -testid t02 -bg cicev6.0.0bs

cd base_suite.t02

# Once all jobs finish, concatenate all output

./results.csh  # All tests results will be stored in results.log

# To plot a timeseries of "total ice extent", "total ice area", and "total ice volume"

./timeseries.csh <directory>

ls \*.png

***********************************************
To run a test suite to compare to baseline data
***********************************************

./cice.setup -m wolf -ts base_suite -testid t03 -bc cicev6.0.0bs

cd base_suite.t03

# Once all jobs finish, concatenate all output

./results.csh  # All tests results will be stored in results.log

# To plot a timeseries of "total ice extent", "total ice area", and "total ice volume"

./timeseries.csh <directory>

ls \*.png

**************************
To compare to another test
**************************
`First:`

./cice.setup -m wolf -t smoke -testid t01 -p 8x2

cd wolf_smoke_gx3_8x2.t01

./cice.build

./cice.submit

# After job finishes, check output

cat test_output

`Then, do the comparison:` 

./cice.setup -m wolf -t smoke -testid t01 -td smoke_gx3_8x2 -s thread -p 4x1

cd wolf_smoke_gx3_4x1_thread.t01

./cice.build

./cice.submit

# After job finishes, check output

cat test_output

******************
Additional Details
******************

- In general, the baseline generation, baseline compare, and test diff are independent.
- Use the '-bd' flag to specify the location where you want the baseline dataset
    to be written.  Without specifying '-bd', the baseline dataset will be written
    to the default baseline directory found in the env.<machine> file (ICE_MACHINE_BASELINE).
- If '-bd' is not passed, the scripts will look for baseline datasets in the default 
    baseline directory found in the env.<machine> file (ICE_MACHINE_BASELINE).
    If the '-bd' option is passed, the scripts will look for baseline datasets in the
    location passed to the -bd argument.
- To generate a baseline dataset for a specific version (for regression testing),
    use '-bg <version_name>'.  The scripts will then place the baseline dataset
    in $ICE_MACHINE_BASELINE/<version_name>/
- The '-testid' flag allows users to specify a testing id that will be added to the
    end of the case directory.  For example, "./cice.setup -m wolf -t smoke -testid t12 -p 4x1"
    creates the directory wolf_smoke_gx3_4x1.t12.  This flag is REQUIRED if using -t or -ts.

.. _compliance:

~~~~~~~~~~~~~~~~~~~~
Code Compliance Test
~~~~~~~~~~~~~~~~~~~~

A core tenet of CICE dycore and Icepack innovations is that they must not change 
the physics and biogeochemistry of existing model configurations, notwithstanding 
obsolete model components. Therefore, alterations to existing CICE Consortium code
must only fix demonstrable numerical or scientific inaccuracies or bugs, or be 
necessary to introduce new science into the code.  New physics and biogeochemistry 
introduced into the model must not change model answers when switched off, and in 
that case CICEcore and Icepack must reproduce answers bit-for-bit as compared to 
previous simulations with the same namelist configurations. This bit-for-bit 
requirement is common in Earth System Modeling projects, but often cannot be achieved 
in practice because model additions may require changes to existing code.  In this 
circumstance, bit-for-bit reproducibility using one compiler may not be unachievable 
on a different computing platform with a different compiler.  Therefore, tools for 
scientific testing of CICE code changes have been developed to accompany bit-for-bit 
testing. These tools exploit the statistical properties of simulated sea ice thickness 
to confirm or deny the null hypothesis, which is that new additions to the CICE dycore 
and Icepack have not significantly altered simulated ice volume using previous model 
configurations.  Here we describe the CICE testing tools, which are applies to output 
from five-year gx-1 simulations that use the standard CICE atmospheric forcing. 
A scientific justification of the testing is provided in
:cite:`Hunke2018`.

.. _paired:

*******************************
Two-Stage Paired Thickness Test
*******************************

The first quality check aims to confirm the null hypotheses
:math:`H_0\!:\!\mu_d{=}0` at every model grid point, given the mean
thickness difference :math:`\mu_d` between paired CICE simulations
‘:math:`a`’ and ‘:math:`b`’ that should be identical. :math:`\mu_d` is
approximated as
:math:`\bar{h}_{d}=\tfrac{1}{n}\sum_{i=1}^n (h_{ai}{-}h_{bi})` for
:math:`n` paired samples of ice thickness :math:`h_{ai}` and
:math:`h_{bi}` in each grid cell of the gx-1 mesh. Following
:cite:`Wilks2006`, the associated :math:`t`-statistic
expects a zero mean, and is therefore

.. math::
   t=\frac{\bar{h}_{d}}{\sigma_d/\sqrt{n_{eff}}}
   :label: t-distribution

given variance
:math:`\sigma_d^{\;2}=\frac{1}{n-1}\sum_{i=1}^{n}(h_{di}-\bar{h}_d)^2`
of :math:`h_{di}{=}(h_{ai}{-}h_{bi})` and effective sample size

.. math::
   n_{eff}{=}n\frac{({1-r_1})}{({1+r_1})}
   :label: neff

for lag-1 autocorrelation:

.. math::
   r_1=\frac{\sum\limits_{i=1}^{n-1}\big[(h_{di}-\bar{h}_{d1:n-1})(h_{di+1}-\bar{h}_{d2:n})\big]}{\sqrt{\sum\limits_{i=1}^{n-1} (h_{di}-\bar{h}_{d1:n-1})^2 \sum\limits_{i=2}^{n} (h_{di}-\bar{h}_{d2:n})^2 }}.
   :label: r1

Here, :math:`\bar{h}_{d1:n-1}` is the mean of all samples except the
last, and :math:`\bar{h}_{d2:n}` is the mean of samples except the
first, and both differ from the overall mean :math:`\bar{h}_d` in
equations (:eq:`t-distribution`). That is:

.. math::
   \bar{h}_{d1:n-1}=\frac{1}{n{-}1} \sum \limits_{i=1}^{n-1} h_{di},\quad 
   \bar{h}_{d2:n}=\frac{1}{n{-}1} \sum \limits_{i=2}^{n} h_{di},\quad
   \bar{h}_d=\frac{1}{n} \sum \limits_{i=1}^{n} {h}_{di}
   :label: short-means

Following :cite:`Zwiers1995`, the effective sample size is
limited to :math:`n_{eff}\in[2,n]`. This definition of :math:`n_{eff}`
assumes ice thickness evolves as an AR(1) process
:cite:`VonStorch1999`, which can be justified by analyzing
the spectral density of daily samples of ice thickness from 5-year
records in CICE Consortium member models :cite:`Hunke2018`.
The AR(1) approximation is inadmissible for paired velocity samples,
because ice drift possesses periodicity from inertia and tides
:cite:`Hibler2006,Lepparanta2012,Roberts2015`. Conversely,
tests of paired ice concentration samples may be less sensitive to ice
drift than ice thickness. In short, ice thickness is the best variable
for CICE Consortium quality control (QC), and for the test of the mean
in particular.

Care is required in analyzing mean sea ice thickness changes using
(:eq:`t-distribution`) with
:math:`N{=}n_{eff}{-}1` degrees of freedom.
:cite:`Zwiers1995` demonstrate that the :math:`t`-test in
(:eq:`t-distribution`) becomes conservative when
:math:`n_{eff} < 30`, meaning that :math:`H_0` may be erroneously
confirmed for highly auto-correlated series. Strong autocorrelation
frequently occurs in modeled sea ice thickness, and :math:`r_1>0.99` is
possible in parts of the gx-1 domain for the five-year QC simulations.
In the event that :math:`H_0` is confirmed but :math:`2\leq n_{eff}<30`,
the :math:`t`-test progresses to the ‘Table Lookup Test’ of
:cite:`Zwiers1995`, to check that the first-stage test
using (:eq:`t-distribution`) was not
conservative. The Table Lookup Test chooses critical :math:`t` values
:math:`|t|<t_{crit}({1{-}\alpha/2},N)` at the :math:`\alpha`
significance level based on :math:`r_1`. It uses the conventional
:math:`t={\bar{h}_{d} \sqrt{n}}/{\sigma_d}` statistic with degrees of
freedom :math:`N{=}n{-}1`, but with :math:`t_{crit}` values generated
using the Monte Carlo technique described in
:cite:`Zwiers1995`, and summarized in :ref:`Table-Lookup` for 5-year QC
simulations (:math:`N=1824`) at the two-sided 80% confidence interval
(:math:`\alpha=0.2`). We choose this interval to limit Type II errors,
whereby a QC test erroneously confirms :math:`H_0`.

:ref:`Table-Lookup` : Summary of two-sided :math:`t_{crit}` values for the Table
Lookup Test of :cite:`Zwiers1995` at the 80% confidence
interval generated for :math:`N=1824` degrees of freedom and lag-1
autocorrelation :math:`r_1`.

.. _Table-Lookup:

.. csv-table:: Table 1
   :widths: 10, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5

   :math:`r_1`,-0.05,0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.97,0.99
   :math:`t_{crit}`,1.32,1.32,1.54,2.02,2.29,2.46,3.17,3.99,5.59,8.44,10.85,20.44

| 
| 

.. _quadratic:

*******************************
Quadratic Skill Compliance Test
*******************************

In addition to the two-stage test of mean sea ice thickness, we also
check that paired simulations are highly correlated and have similar
variance using a skill metric adapted from
:cite:`Taylor2001`. A general skill score applicable to
Taylor diagrams takes the form

.. math::
   S_m=\frac{4(1+R)^m}{({\hat{\sigma}_{f}+1/{\hat{\sigma}_{f}}})^2 (1+R_0)^m}
   :label: taylor-skill

where :math:`m=1` for variance-weighted skill, and :math:`m=4` for
correlation-weighted performance, as given in equations (4) and (5) of
:cite:`Taylor2001`, respectively. We choose :math:`m=2` to
balance the importance of variance and correlation reproduction in QC
tests, where :math:`\hat{\sigma}_{f}={\sigma_{b}}/{\sigma_{a}}` is the ratio
of the standard deviations of simulations ‘:math:`b`’ and ‘:math:`a`’,
respectively, and simulation ‘:math:`a`’ is the control. :math:`R_0` is
the maximum possible correlation between two series for correlation
coefficient :math:`R` calculated between respective thickness pairs
:math:`h_{a}` and :math:`h_{b}`. Bit-for-bit reproduction of previous
CICE simulations means that perfect correlation is possible, and so
:math:`R_0=1`, giving the quadratic skill of run ‘:math:`b`’ relative to
run ‘:math:`a`’:

.. math::
   S=\bigg[ \frac{(1+R) (\sigma_a \sigma_b)}{({\sigma_a}^2 + {\sigma_b}^2)} \bigg]^2
   :label: quadratic-skill

This provides a skill score between 0 and 1. We apply this :math:`S`
metric separately to the northern and southern hemispheres of the gx-1
grid by area-weighting the daily thickness samples discussed in the
Two-Stage Paired Thickness QC Test. The hemispheric mean thickness over
a 5-year simulation for run ‘:math:`a`’ is:

.. math::
   \bar{h}_{a}=\frac{1}{n} \sum_{i=1}^{n} \sum_{j=1}^{J} \ W_{j} \; h_{{a}_{i,j}}
   :label: h-bar

at time sample :math:`i` and grid point index :math:`j`, with an
equivalent equation for simulation ‘:math:`b`’. :math:`n` is the total
number of time samples (nominally :math:`n=1825`) and :math:`J` is the
total number of grid points on the gx-1 grid. :math:`W_j` is the weight
attributed to each grid point according to its area :math:`A_{j}`, given
as

.. math::
   W_{j}=\frac{ A_{j} }{\sum_{j=1}^{J} A_{j}}
   :label: area-weight

for all grid points within each hemisphere with one or more non-zero
thicknesses in one or both sets of samples :math:`h_{{a}_{i,j}}` or
:math:`h_{{b}_{i,j}}`. The area-weighted variance for simulation
‘:math:`a`’ is:

.. math::
   \sigma_a^{\;2}=\frac{\hat{J}}{(n\,\hat{J}-1)} \sum_{i=1}^{n} \sum_{j=1}^{J}  W_{j} \, (h_{{a}_{i,j}}-\bar{h}_{a})^2
   :label: weighted-deviation

where :math:`\hat{J}` is the number of non-zero :math:`W_j` weights,
and :math:`\sigma_b` is calculated equivalently for run ‘:math:`b`’. In
this context, :math:`R` becomes a weighted correlation coefficient,
calculated as

.. math::
   R=\frac{\textrm{cov}(h_{a},h_{b})}{\sigma_a \; \sigma_b}
   :label: R

given the weighted covariance

.. math::
   \textrm{cov}(h_{a},h_{b})=\frac{\hat{J}}{(n\,\hat{J}-1)} \sum_{i=1}^{n} \sum_{j=1}^{J}  W_{j} \, (h_{{a}_{i,j}}-\bar{h}_{a}) (h_{{b}_{i,j}}-\bar{h}_{b}).
   :label: weighted-covariance

Using equations (:eq:`quadratic-skill`)
to (:eq:`weighted-covariance`), the skill
score :math:`S` is calculated separately for the northern and southern
hemispheres, and must exceed a critical value nominally set to
:math:`S_{crit}=0.99` to pass the test. Practical illustrations of this
test and the Two-Stage test described in the previous section are
provided in :cite:`Hunke2018`.

***************************
Practical Testing Procedure
***************************

The CICE code compliance test is performed by running a python script (cice.t-test.py).
In order to run the script, the following requirements must be met:

* Python v2.7 or later
* netCDF Python package
* numpy Python package
* matplotlib Python package (optional)
* basemap Python package (optional)

In order to generate the files necessary for the compliance test, test cases should be
created with the ``qc`` option (i.e., ``-s qc``) when running cice.setup.  This 
option results in daily, non-averaged history files being written for a 5 year simulation.

To install the necessary Python packages, the ``pip`` Python utility can be used.

.. code-block:: bash

  pip install --user netCDF4
  pip install --user numpy
  pip install --user matplotlib

To run the compliance test:

.. code-block:: bash

  cp configuration/scripts/tests/QC/cice.t-test.py .
  ./cice.t-test.py /path/to/baseline/history /path/to/test/history

The script will produce output similar to:

  |  \INFO:__main__:Number of files: 1825
  |  \INFO:__main__:Two-Stage Test Passed
  |  \INFO:__main__:Quadratic Skill Test Passed for Northern Hemisphere
  |  \INFO:__main__:Quadratic Skill Test Passed for Southern Hemisphere
  |  \INFO:__main__:
  |  \INFO:__main__:Quality Control Test PASSED

Additionally, the exit code from the test (``echo $?``) will be 0 if the test passed,
and 1 if the test failed.

Implementation notes: 1) Provide a pass/fail on each of the confidence
intervals, 2) Facilitate output of a bitmap for each test so that
locations of failures can be identified.

~~~~~~~~~~~~~~~~~~~~
CICE Test Reporting
~~~~~~~~~~~~~~~~~~~~

The CICE testing scripts have the capability of posting the test results 
to an online dashboard, located `on CDash <http://my.cdash.org/index.php?project=myCICE>`_.  
There are 2 options for posting CICE results to CDash: 1) The automated 
script, 2) The manual method.

*****************
Automatic Script
*****************

To automatically run the CICE tests, and post the results to the CICE Cdash dashboard,
users need to copy and run the ``run.suite`` script:

.. code-block:: bash

  cp configuration/scripts/run.suite .
  ./run.suite -m <machine> -testid <test_id> -bc <baseline_to_compare> -bg <baseline_to_generate>

The run.suite script does the following:

- Creates a fresh clone of the CICE-Consortium repository
- ``cd`` to cloned repo
- run ``cice.setup`` to generate the base_suite directories.  The output 
  is piped to ``log.suite``
- Running ``cice.setup`` submits each individual job to the queue.  
- ``run.suite`` monitors the queue manager to determine when all jobs have 
  finished (pings the queue manager once every 5 minutes).
- Once all jobs complete, cd to base_suite directory and run ``./results.csh``
- Run ``./run_ctest.csh`` in order to post the test results to the CDash dashboard

*****************
Manual Method
*****************

To manually run the CICE tests and post the results to the CICE CDash dashboard,
users essentially just need to perform all steps available in run.suite, detailed below:

- Pass the ``-report`` flag to cice.setup when running the ``base_suite`` test suite.
  The ``-report`` flag copies the required CTest / CDash scripts to the suite 
  directory.
- ``cice.setup`` compiles the CICE code, and submits all of the jobs to the 
  queue manager.  
- After every job has been submitted and completed, ``cd`` to the suite directory.
- Parse the results, by running ``./results.csh``.
- Run the CTest / CDash script ``./run_ctest.csh``.

If the ``run_ctest.csh`` script is unable to post the testing results to the CDash
server, a message will be printed to the screen detailing instructions on how to attempt
to post the results from another server.  If ``run_ctest.csh`` fails to submit the results,
it will generate a tarball ``cice_ctest.tgz`` that contains the necessary files for 
submission.  Copy this file to another server (CMake version 2.8+ required), extract the 
archive, and run ``./run_ctest.csh -submit``.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~
End-To-End Testing Procedure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example of a step-by-step procedure for testing a code change that results
in non-bit-for-bit results:

.. code-block:: bash

  # Create a baseline dataset (only necessary if no baseline exists on the system)
  ./cice.setup -m onyx -ts base_suite -testid base0 -bg cicev6.0.0 -a <account_number>

  # Check out the updated code, or clone from a pull request

  # Run the test with the new code
  ./cice.setup -m onyx -ts base_suite -testid test0 -bc cicev6.0.0 -a <account_number>

  # Check the results
  cd base_suite.test0
  ./results.csh

  #### If the BFB tests fail, perform the compliance testing ####
  # Create a QC baseline
  ./cice.setup -m onyx -t smoke -g gx1 -p 44x1 -testid qc_base -s qc,medium -a <account_number>
  cd onyx_smoke_gx1_44x1_medium_qc.qc_base
  ./cice.build
  ./cice.submit

  # Check out the updated code or clone from a pull request

  # Create the t-test testing data
  ./cice.setup -m onyx -t smoke -g gx1 -p 44x1 -testid qc_test -s qc,medium -a <account_number>
  cd onyx_smoke_gx1_44x1_medium_qc.qc_test
  ./cice.build
  ./cice.submit

  # Wait for runs to finish
  
  # Perform the QC test
  cp configuration/scripts/tests/QC/cice.t-test.py
  ./cice.t-test.py /p/work/turner/CICE_RUNS/onyx_smoke_gx1_44x1_medium_qc.qc_base \
                   /p/work/turner/CICE_RUNS/onyx_smoke_gx1_44x1_medium_qc.qc_test

  # Example output:
  INFO:__main__:Number of files: 1825
  INFO:__main__:Two-Stage Test Passed
  INFO:__main__:Quadratic Skill Test Passed for Northern Hemisphere
  INFO:__main__:Quadratic Skill Test Passed for Southern Hemisphere
  INFO:__main__:
  INFO:__main__:Quality Control Test PASSED

.. _tabnamelist:

-------------------------
Table of namelist options
-------------------------

.. _tab-namelist:

.. csv-table:: Table 8
   :header: "variable", "options/format", "description", "recommended value"
   :widths: 15, 15, 30, 15 

   "*setup_nml*", " ", " ", " "
   "", "", "*Time, Diagnostics*", ""
   "``days_per_year``", "``360`` or ``365``", "number of days in a model year", "365"
   "``use_leap_years``", "true/false", "if true, include leap days", ""
   "``year_init``", "yyyy", "the initial year, if not using restart", ""
   "``istep0``", "integer", "initial time step number", "0"
   "``dt``", "seconds", "thermodynamics time step length", "3600."
   "``npt``", "integer", "total number of time steps to take", ""
   "``ndtd``", "integer", "number of dynamics/advection/ridging/steps per thermo timestep", "1"
   "", "", "*Initialization/Restarting*", ""
   "``runtype``", "``initial``", "start from ``ice_ic``", ""
   "", "``continue``", "restart using ``pointer_file``", ""
   "``ice_ic``", "``default``", "latitude and sst dependent", "default"
   "", "``none``", "no ice", ""
   "", "path/file", "restart file name", ""
   "``restart``", "true/false", "initialize using restart file", "``.true.``"
   "``use_restart_time``", "true/false", "set initial date using restart file", "``.true.``"
   "``restart_format``", "nc", "read/write  restart files (use with PIO)", ""
   "", "bin", "read/write binary restart files", ""
   "``lcdf64``", "true/false", "if true, use 64-bit  format", ""
   "``restart_dir``", "path/", "path to restart directory", ""
   "``restart_ext``", "true/false", "read/write halo cells in restart files", ""
   "``restart_file``", "filename prefix", "output file for restart dump", "‘iced’"
   "``pointer_file``", "pointer filename", "contains restart filename", ""
   "``dumpfreq``", "``y``", "write restart every ``dumpfreq_n`` years", "y"
   "", "``m``", "write restart every ``dumpfreq_n`` months", ""
   "", "``d``", "write restart every ``dumpfreq_n`` days", ""
   "``dumpfreq_n``", "integer", "frequency restart data is written", "1"
   "``dump_last``", "true/false", "if true, write restart on last time step of simulation", ""
   "", "", "*Model Output*", ""
   "``bfbflag``", "true/false", "for bit-for-bit diagnostic output", ""
   "``diagfreq``", "integer", "frequency of diagnostic output in ``dt``", "24"
   "", "*e.g.*, 10", "once every 10 time steps", ""
   "``diag_type``", "``stdout``", "write diagnostic output to stdout", ""
   "", "``file``", "write diagnostic output to file", ""
   "``diag_file``", "filename", "diagnostic output file (script may reset)", ""
   "``print_global``", "true/false", "print diagnostic data, global sums", "``.false.``"
   "``print_points``", "true/false", "print diagnostic data for two grid points", "``.false.``"
   "``latpnt``", "real", "latitude of (2) diagnostic points", "" 
   "``lonpnt``", "real", "longitude of (2) diagnostic points", ""
   "``dbug``", "true/false", "if true, write extra diagnostics", "``.false.``"
   "``histfreq``", "string array", "defines output frequencies", ""
   "", "``y``", "write history every ``histfreq_n`` years", ""
   "", "``m``", "write history every ``histfreq_n`` months", ""
   "", "``d``", "write history every ``histfreq_n`` days", ""
   "", "``h``", "write history every ``histfreq_n`` hours", ""
   "", "``1``", "write history every time step", ""
   "", "``x``", "unused frequency stream (not written)", ""
   "``histfreq_n``", "integer array", "frequency history output is written", ""
   "", "0", "do not write to history", ""
   "``hist_avg``", "true", "write time-averaged data", "``.true.``"
   "", "false", "write snapshots of data", ""
   "``history_dir``", "path/", "path to history output directory", ""
   "``history_file``", "filename prefix", "output file for history", "‘iceh’"
   "``write_ic``", "true/false", "write initial condition", ""
   "``incond_dir``", "path/", "path to initial condition directory", ""
   "``incond_file``", "filename prefix", "output file for initial condition", "‘iceh’"
   "``runid``", "string", "label for run (currently CESM only)", ""
   "", "", "", ""
   "*grid_nml*", "", "", ""
   "", "", "*Grid*", ""
   "``grid_format``", "``nc``", "read  grid and kmt files", "‘bin’"
   "", "``bin``", "read direct access, binary file", ""
   "``grid_type``", "``rectangular``", "defined in *rectgrid*", ""
   "", "``displaced_pole``", "read from file in *popgrid*", ""
   "", "``tripole``", "read from file in *popgrid*", ""
   "", "``regional``", "read from file in *popgrid*", ""
   "``grid_file``", "filename", "name of grid file to be read", "‘grid’"
   "``kmt_file``", "filename", "name of land mask file to be read", "‘kmt’"
   "``gridcpl_file``", "filename", "input file for coupling grid info", ""
   "``kcatbound``", "``0``", "original category boundary formula", "0"
   "", "``1``", "new formula with round numbers", ""
   "", "``2``", "WMO standard categories", ""
   "", "``-1``", "one category", ""
   "", "", "", ""
   "*domain_nml*", "", "", ""
   "", "", "*Domain*", ""
   "``nprocs``", "integer", "number of processors to use", ""
   "``processor_shape``", "``slenderX1``", "1 processor in the y direction (tall, thin)", ""
   "", "``slenderX2``", "2 processors in the y direction (thin)", ""
   "", "``square-ice``", "more processors in x than y, :math:`\sim` square", ""
   "", "``square-pop``", "more processors in y than x, :math:`\sim` square", ""
   "``distribution_type``", "``cartesian``", "distribute blocks in 2D Cartesian array", ""
   "", "``roundrobin``", "1 block per proc until blocks are used", ""
   "", "``sectcart``", "blocks distributed to domain quadrants", ""
   "", "``sectrobin``", "several blocks per proc until used", ""
   "", "``rake``", "redistribute blocks among neighbors", ""
   "", "``spacecurve``", "distribute blocks via space-filling curves", ""
   "``distribution_weight``", "``block``", "full block size sets ``work_per_block``", ""
   "", "``latitude``", "latitude/ocean sets ``work_per_block``", ""
   "``ew_boundary_type``", "``cyclic``", "periodic boundary conditions in x-direction", ""
   "", "``open``", "Dirichlet boundary conditions in x", ""
   "``ns_boundary_type``", "``cyclic``", "periodic boundary conditions in y-direction", ""
   "", "``open``", "Dirichlet boundary conditions in y", ""
   "", "``tripole``", "U-fold tripole boundary conditions in y", ""
   "", "``tripoleT``", "T-fold tripole boundary conditions in y", ""
   "``maskhalo_dyn``", "true/false", "mask unused halo cells for dynamics", ""
   "``maskhalo_remap``", "true/false", "mask unused halo cells for transport", ""
   "``maskhalo_bound``", "true/false", "mask unused halo cells for boundary updates", ""
   "", "", "", ""
   "*tracer_nml*", "", "", ""
   "", "", "*Tracers*", ""
   "``tr_iage``", "true/false", "ice age", ""
   "``restart_age``", "true/false", "restart tracer values from file", ""
   "``tr_FY``", "true/false", "first-year ice area", ""
   "``restart_FY``", "true/false", "restart tracer values from file", ""
   "``tr_lvl``", "true/false", "level ice area and volume", ""
   "``restart_lvl``", "true/false", "restart tracer values from file", ""
   "``tr_pond_cesm``", "true/false", "CESM melt ponds", ""
   "``restart_pond_cesm``", "true/false", "restart tracer values from file", ""
   "``tr_pond_topo``", "true/false", "topo melt ponds", ""
   "``restart_pond_topo``", "true/false", "restart tracer values from file", ""
   "``tr_pond_lvl``", "true/false", "level-ice melt ponds", ""
   "``restart_pond_lvl``", "true/false", "restart tracer values from file", ""
   "``tr_aero``", "true/false", "aerosols", ""
   "``restart_aero``", "true/false", "restart tracer values from file", ""
   "*thermo_nml*", "", "", ""
   "", "", "*Thermodynamics*", ""
   "``kitd``", "``0``", "delta function ITD approximation", "1"
   "", "``1``", "linear remapping ITD approximation", ""
   "``ktherm``", "``0``", "zero-layer thermodynamic model", ""
   "", "``1``", "Bitz and Lipscomb thermodynamic model", ""
   "", "``2``", "mushy-layer thermodynamic model", ""
   "``conduct``", "``MU71``", "conductivity :cite:`MU71`", ""
   "", "``bubbly``", "conductivity :cite:`PETB07`", ""
   "``a_rapid_mode``", "real", "brine channel diameter", "0.5x10 :math:`^{-3}` m"
   "``Rac_rapid_mode``", "real", "critical Rayleigh number", "10"
   "``aspect_rapid_mode``", "real", "brine convection aspect ratio", "1"
   "``dSdt_slow_mode``", "real", "drainage strength parameter", "-1.5x10 :math:`^{-7}` m/s/K"
   "``phi_c_slow_mode``", ":math:`0<\phi_c < 1`", "critical liquid fraction", "0.05"
   "``phi_i_mushy``", ":math:`0<\phi_i < 1`", "solid fraction at lower boundary", "0.85"
   "", "", "", ""
   "*dynamics_nml*", "", "", ""
   "", "", "*Dynamics*", ""
   "``kdyn``", "``0``", "dynamics OFF", "1"
   "", "``1``", "EVP dynamics", ""
   "", "``2``", "EAP dynamics", ""
   "``revised_evp``", "true/false", "use revised EVP formulation", ""
   "``ndte``", "integer", "number of EVP subcycles", "120"
   "``advection``", "``remap``", "linear remapping advection", "‘remap’"
   "", "``upwind``", "donor cell advection", ""
   "``kstrength``", "``0``", "ice strength formulation :cite:`Hibler79`", "1"
   "", "``1``", "ice strength formulation :cite:`Rothrock75`", ""
   "``krdg_partic``", "``0``", "old ridging participation function", "1"
   "", "``1``", "new ridging participation function", ""
   "``krdg_redist``", "``0``", "old ridging redistribution function", "1"
   "", "``1``", "new ridging redistribution function", ""
   "``mu_rdg``", "real", "e-folding scale of ridged ice", ""
   "``Cf``", "real", "ratio of ridging work to PE change in ridging", "17."
   "", "", "", ""
   "*shortwave_nml*", "", "", ""
   "", "", "*Shortwave*", ""
   "``shortwave``", "``default``", "NCAR CCSM3 distribution method", ""
   "", "``dEdd``", "Delta-Eddington method", ""
   "``albedo_type``", "``default``", "NCAR CCSM3 albedos", "‘default’"
   "", "``constant``", "four constant albedos", ""
   "``albicev``", ":math:`0<\alpha <1`", "visible ice albedo for thicker ice", ""
   "``albicei``", ":math:`0<\alpha <1`", "near infrared ice albedo for thicker ice", ""
   "``albsnowv``", ":math:`0<\alpha <1`", "visible, cold snow albedo", ""
   "``albsnowi``", ":math:`0<\alpha <1`", "near infrared, cold snow albedo", ""
   "``ahmax``", "real", "albedo is constant above this thickness", "0.3 m"
   "``R_ice``", "real", "tuning parameter for sea ice albedo from Delta-Eddington shortwave", ""
   "``R_pnd``", "real", "... for ponded sea ice albedo …", ""
   "``R_snw``", "real", "... for snow (broadband albedo) …", ""
   "``dT_mlt``", "real", ":math:`\Delta` temperature per :math:`\Delta` snow grain radius", ""
   "``rsnw_mlt``", "real", "maximum melting snow grain radius", ""
   "``kalg``", "real", "absorption coefficient for algae", ""
   "", "", "", ""
   "*ponds_nml*", "", "", ""
   "", "", "*Melt Ponds*", ""
   "``hp1``", "real", "critical ice lid thickness for topo ponds", "0.01 m"
   "``hs0``", "real", "snow depth of transition to bare sea ice", "0.03 m"
   "``hs1``", "real", "snow depth of transition to pond ice", "0.03 m"
   "``dpscale``", "real", "time scale for flushing in permeable ice", ":math:`1\times 10^{-3}`"
   "``frzpnd``", "``hlid``", "Stefan refreezing with pond ice thickness", "‘hlid’"
   "", "``cesm``", "CESM refreezing empirical formula", ""
   "``rfracmin``", ":math:`0 \le r_{min} \le 1`", "minimum melt water added to ponds", "0.15"
   "``rfracmax``", ":math:`0 \le r_{max} \le 1`", "maximum melt water added to ponds", "1.0"
   "``pndaspect``", "real", "aspect ratio of pond changes (depth:area)", "0.8"
   "", "", "", ""
   "*zbgc_nml*", "", "", ""
   "", "", "*Biogeochemistry*", ""
   "``tr_brine``", "true/false", "brine height tracer", ""
   "``tr_zaero``", "true/false", "vertical aerosol tracers", ""
   "``modal_aero``", "true/false", "modal aersols", ""
   "``restore_bgc``", "true/false", "restore bgc to data", ""
   "``solve_zsal`", "true/false", "update salinity tracer profile", ""
   "``bgc_data_dir``", "path/", "data directory for bgc", ""
   "``skl_bgc``", "true/false", "biogeochemistry", ""
   "``sil_data_type``", "``default``", "default forcing value for silicate", ""
   "", "``clim``", "silicate forcing from ocean climatology :cite:`GLBA06`", ""
   "``nit_data_type``", "``default``", "default forcing value for nitrate", ""
   "", "``clim``", "nitrate forcing from ocean climatology :cite:`GLBA06`", ""
   "", "``sss``", "nitrate forcing equals salinity", ""
   "``fe_data_type``", "``default``", "default forcing value for iron", ""
   "", "``clim``", "iron forcing from ocean climatology", ""
   "``bgc_flux_type``", "``Jin2006``", "ice–ocean flux velocity of :cite:`JDWSTWLG06`", ""
   "", "``constant``", "constant ice–ocean flux velocity", ""
   "``restart_bgc``", "true/false", "restart tracer values from file", ""
   "``tr_bgc_C_sk``", "true/false", "algal carbon tracer", ""
   "``tr_bgc_chl_sk``", "true/false", "algal chlorophyll tracer", ""
   "``tr_bgc_Am_sk``", "true/false", "ammonium tracer", ""
   "``tr_bgc_Sil_sk``", "true/false", "silicate tracer", ""
   "``tr_bgc_DMSPp_sk``", "true/false", "particulate DMSP tracer", ""
   "``tr_bgc_DMSPd_sk``", "true/false", "dissolved DMSP tracer", ""
   "``tr_bgc_DMS_sk``", "true/false", "DMS tracer", ""
   "``phi_snow``", "real", "snow porosity for brine height tracer", ""
   "", "", "", ""
   "*forcing_nml*", "", "", ""
   "", "", "*Forcing*", ""
   "``formdrag``", "true/false", "calculate form drag", ""
   "``atmbndy``", "``default``", "stability-based boundary layer", "‘default’"
   "", "``constant``", "bulk transfer coefficients", ""
   "``fyear_init``", "yyyy", "first year of atmospheric forcing data", ""
   "``ycycle``", "integer", "number of years in forcing data cycle", ""
   "``atm_data_format``", "``nc``", "read  atmo forcing files", ""
   "", "``bin``", "read direct access, binary files", ""
   "``atm_data_type``", "``default``", "constant values defined in the code", ""
   "", "``LYq``", "AOMIP/Large-Yeager forcing data", ""
   "", "``monthly``", "monthly forcing data", ""
   "", "``ncar``", "NCAR bulk forcing data", ""
   "", "``oned``", "column forcing data", ""
   "``atm_data_dir``", "path/", "path to atmospheric forcing data directory", ""
   "``calc_strair``", "true", "calculate wind stress and speed", ""
   "", "false", "read wind stress and speed from files", ""
   "``highfreq``", "true/false", "high-frequency atmo coupling", ""
   "``natmiter``", "integer", "number of atmo boundary layer iterations", ""
   "``calc_Tsfc``", "true/false", "calculate surface temperature", "``.true.``"
   "``precip_units``", "``mks``", "liquid precipitation data units", ""
   "", "``mm_per_month``", "", ""
   "", "``mm_per_sec``", "(same as MKS units)", ""
   "``tfrz_option``", "``minus1p8``", "constant ocean freezing temperature (:math:`-1.8^{\circ} C`)", ""
   "", "``linear_salt``", "linear function of salinity (ktherm=1)", ""
   "", "``mushy_layer``", "matches mushy-layer thermo (ktherm=2)", ""
   "``ustar_min``", "real", "minimum value of ocean friction velocity", "0.0005 m/s"
   "``fbot_xfer_type``", "``constant``", "constant ocean heat transfer coefficient", ""
   "", "``Cdn_ocn``", "variable ocean heat transfer coefficient", ""
   "``update_ocn_f``", "true", "include frazil water/salt fluxes in ocn fluxes", ""
   "", "false", "do not include (when coupling with POP)", ""
   "``l_mpond_fresh``", "true", "retain (topo) pond water until ponds drain", ""
   "", "false", "release (topo) pond water immediately to ocean", ""
   "``oceanmixed_ice``", "true/false", "active ocean mixed layer calculation", "``.true.`` (if uncoupled)"
   "``ocn_data_format``", "``nc``", "read  ocean forcing files", ""
   "", "``bin``", "read direct access, binary files", ""
   "``sss_data_type``", "``default``", "constant values defined in the code", ""
   "", "``clim``", "climatological data", ""
   "", "``near``", "POP ocean forcing data", ""
   "``sst_data_type``", "``default``", "constant values defined in the code", ""
   "", "``clim``", "climatological data", ""
   "", "``ncar``", "POP ocean forcing data", ""
   "``ocn_data_dir``", "path/", "path to oceanic forcing data directory", ""
   "``oceanmixed_file``", "filename", "data file containing ocean forcing data", ""
   "``restore_sst``", "true/false", "restore sst to data", ""
   "``trestore``", "integer", "sst restoring time scale (days)", ""
   "``restore_ice``", "true/false", "restore ice state along lateral boundaries", ""
   "", "", "", ""
   "*icefields_tracer_nml*", "", "", ""
   "", "", "*History Fields*", ""
   "``f_<var>``", "string", "frequency units for writing ``<var>`` to history", ""
   "", "``y``", "write history every ``histfreq_n`` years", ""
   "", "``m``", "write history every ``histfreq_n`` months", ""
   "", "``d``", "write history every ``histfreq_n`` days", ""
   "", "``h``", "write history every ``histfreq_n`` hours", ""
   "", "``1``", "write history every time step", ""
   "", "``x``", "do not write ``<var>`` to history", ""
   "", "``md``", "*e.g.,* write both monthly and daily files", ""
   "``f_<var>_ai``", "", "grid cell average of ``<var>`` (:math:`\times a_i`)", ""

