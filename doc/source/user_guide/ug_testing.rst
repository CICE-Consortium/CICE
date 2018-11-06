:tocdepth: 3

.. _testing:

Testing CICE
================

This section documents primarily how to use the CICE scripts to carry 
out CICE testing.  Exactly what to test is a separate question and
depends on the kinds of code changes being made.  Prior to merging
changes to the CICE Consortium master, changes will be reviewed and
developers will need to provide a summary of the tests carried out.

There is a base suite of tests provided by default with CICE and this
may be a good starting point for testing.

The testing scripts support several features
 - Ability to test individual (via ``--test``)or multiple tests (via ``--suite``)
   using an input file to define the suite
 - Ability to use test suites defined in the package or test suites defined by the user
 - Ability to store test results for regresssion testing (``--bgen``)
 - Ability to compare results to prior baselines to verify bit-for-bit (``--bcmp``)
 - Ability to define where baseline tests are stored (``--bdir``)
 - Ability to compare tests against each other (``--diff``)
 - Ability to set account number (``--acct``), which is otherwise not set and may result in tests not being submitted

.. _indtests:

Individual Tests
----------------

The CICE scripts support both setup of individual tests as well as test suites.  Individual
tests are run from the command line::

  ./cice.setup --test smoke --mach conrad --env cray --set diag1,debug --testid myid 

Tests are just like cases but have some additional scripting around them.  Individual
tests can be created and manually modified just like cases.
Many of the command line arguments for individual tests
are similar to :ref:`case_options` for ``--case``.  
For individual tests, the following command line options can be set

``--test`` TESTNAME
     specifies the test type.  This is probably either smoke or restart but see `cice.setup --help` for the latest.  This is required instead of ``--case``.

``--testid`` ID
     specifies the testid.  This is required for every use of ``--test`` and ``--suite``.  This is a user defined string that will allow each test to have a unique case and run directory name.  This is also required.

``--mach`` MACHINE (see :ref:`case_options`)

``--env`` ENVIRONMENT1 (see :ref:`case_options`)

``--set`` SET1,SET2,SET3 (see :ref:`case_options`)

``--acct`` ACCOUNT (see :ref:`case_options`)

``--grid`` GRID (see :ref:`case_options`)

``--pes`` MxNxBXxBYxMB (see :ref:`case_options`)

There are several additional options that come with ``--test`` that are not available
with ``--case`` for regression and comparision testing,

``--bdir`` DIR
     specifies the top level location of the baseline results.  This is used in conjuction with ``--bgen`` and ``--bcmp``.  The default is set by ICE_MACHINE_BASELINE in the env.[machine]_[environment] file.

``--bgen`` DIR
     specifies the name of the directory under [bdir] where test results will be stored.  When this flag is set, it automatically creates that directory and stores results from the test under that directory.  If DIR is set to ``default``, then the scripts will automatically generate a directory name based on the CICE hash and the date and time.  This can be useful for tracking the baselines by hash.

``--bcmp`` DIR
     specifies the name of the directory under [bdir] that the current tests will be compared to.  When this flag is set, it automatically invokes regression testing and compares results from the current test to those prior results.  If DIR is set to ``default``, then the script will automatically generate the last directory name in the [bdir] directory.  This can be useful for automated regression testing.

``--diff`` LONG_TESTNAME
     invokes a comparison against another local test.  This allows different tests to be compared to each other for bit-for-bit-ness.  This is different than ``--bcmp``.  ``--bcmp`` is regression testing, comparing identical test results between different model versions.  ``--diff`` allows comparison of two different test cases against each other.  For instance, different block sizes, decompositions, and other model features are expected to produced identical results and ``--diff`` supports that testing.  The restrictions for use of ``--diff`` are that the test has to already be completed and the testid has to match.  The LONG_TESTNAME string should be of format [test]_[grid]_[pes]_[sets].  The [machine], [env], and [testid] will be added to that string to complete the testname being compared.  (See also :ref:`examplediff` #5)

The format of the case directory name for a test will always be 
``[machine]_[env]_[test]_[grid]_[pes]_[sets].[testid]``
The [sets] will always be sorted alphabetically by the script so ``--set debug,diag1`` and
``--set diag1,debug`` produces the same testname and test with _debug_diag1 in that order.

To build and run a test after invoking the ./cice.setup command, the process is the same as for a case.  
cd to the test directory, run the build script, and run the submit script::

 cd [test_case]
 ./cice.build
 ./cice.submit

The test results will be generated in a local file called **test_output**.
To check those results::

 cat test_output

Tests are defined under **configuration/scripts/tests/**.  Some tests currently supported are:

- smoke   - Runs the model for default length.  The length and options can
            be set with the ``--set`` command line option.  The test passes if the
            model completes successfully.
- restart - Runs the model for 10 days, writing a restart file at the end of day 5 and
            again at the end of the run.  Runs the model a second time starting from the
            day 5 restart and writes a restart at then end of day 10 of the model run.
            The test passes if both runs complete and
            if the restart files at the end of day 10 from both runs are bit-for-bit identical.
- decomp   - Runs a set of different decompositions on a given configuration

Please run ``./cice.setup --help`` for the latest information.


Adding a new test
~~~~~~~~~~~~~~~~~~~~~~~~

See :ref:`dev_testing`


.. _examplediff:

Individual Test Examples
~~~~~~~~~~~~~~~~~~~~~~~~

 1) **Basic default single test**
     
    Define the test, mach, env, and testid.
    ::

      ./cice.setup --test smoke --mach wolf --env gnu --testid t00
      cd wolf_gnu_smoke_col_1x1.t00
      ./cice.build
      ./cice.submit
      ./cat test_output

 2) **Simple test with some options**

    Add ``--set``
    ::

      ./cice.setup --test smoke --mach wolf --env gnu --set diag1,debug --testid t00
      cd wolf_gnu_smoke_col_1x1_debug_diag1.t00
      ./cice.build
      ./cice.submit
      ./cat test_output

 3) **Single test, generate a baseline dataset**

    Add ``--bgen`` 
    ::

      ./cice.setup --test smoke --mach wolf -env gnu --bgen cice.v01 --testid t00 --set diag1
      cd wolf_gnu_smoke_col_1x1_diag1.t00
      ./cice.build
      ./cice.submit
      ./cat test_output

 4) **Single test, compare results to a prior baseline**

    Add ``--bcmp``.  For this to work,
    the prior baseline must exist and have the exact same base testname 
    [machine]_[env]_[test]_[grid]_[pes]_[sets] 
    ::

      ./cice.setup --test smoke --mach wolf -env gnu --bcmp cice.v01 --testid t01 --set diag1
      cd wolf_gnu_smoke_col_1x1_diag1.t01
      ./cice.build
      ./cice.submit
      ./cat test_output

 5) **Simple test, generate a baseline dataset and compare to a prior baseline**

    Use ``--bgen`` and ``--bcmp``.  The prior baseline must exist already.
    ::

      ./cice.setup --test smoke --mach wolf -env gnu --bgen cice.v02 --bcmp cice.v01 --testid t02 --set diag1
      cd wolf_gnu_smoke_col_1x1_diag1.t02
      ./cice.build
      ./cice.submit
      ./cat test_output

 6) **Simple test, comparison against another test**

    ``--diff`` provides a way to compare tests with each other.  
    For this to work, the tests have to be run in a specific order and
    the testids need to match.  The test 
    is always compared relative to the current case directory.

    To run the first test,
    ::

      ./cice.setup --test smoke --mach wolf -env gnu --testid tx01 --set debug
      cd wolf_gnu_smoke_col_1x1_debug.tx01
      ./cice.build
      ./cice.submit
      ./cat test_output

    Then to run the second test and compare to the results from the first test
    ::

      ./cice.setup --test smoke --mach wolf -env gnu --testid tx01 --diff smoke_col_1x1_debug
      cd wolf_gnu_smoke_col_1x1.tx01
      ./cice.build
      ./cice.submit
      ./cat test_output

    The scripts will add a [machine]_[environment] to the beginning of the diff 
    argument and the same testid to the end of the diff argument.  Then the runs 
    will be compared for bit-for-bit and a result will be produced in test_output.  

Specific Test Cases
~~~~~~~~~~~~~~~~~~~

In addition to the test implemented in the general testing framework, specific
tests have been developed to validate specific portions of the model.  These
specific tests are detailed in this section.

``box2001``
^^^^^^^^^^^^

The ``box2001`` test case is configured to perform the rectangular-grid box test 
detailed in :cite:`Hunke01`.  It is configured to run a 72-hour simulation with 
thermodynamics disabled in a rectangular domain (80 x 80 grid cells) with a land
boundary around the entire domain.  It includes the following namelist modifications:

- ``dxrect``: ``16.e5`` meters
- ``dyrect``: ``16.e5`` meters
- ``thermo``: ``0`` (disables thermodynamics)
- ``coriolis``: ``zero`` (zero coriolis force)

Ocean stresses are computed as in :cite:`Hunke01` where they are circular and centered 
in the square domain.  The ice distribution is fixed, with a constant 2 meter ice 
thickness and a concentration field that varies linearly in the x-direction from ``0``
to ``1`` and is constant in the y-direction.  No islands are included in this
configuration.  The test is configured to run on a single processor.

To run the test: ``./cice.setup -m <machine> --test smoke -s box2001 --testid <test_id>
--grid gbox80 --acct <queue manager account> -p 1x1``

.. _testsuites:

Test suites
------------

Test suites support running multiple tests specified via
an input file.  When invoking the test suite option (``--suite``) with **cice.setup**,
all tests will be created, built, and submitted automatically under
a local directory called testsuite.[testid] as part of involing the suite.::

  ./cice.setup --suite base_suite --mach wolf --env gnu --testid myid

Like an individual test, the ``--testid`` option must be specified and can be any 
string.  Once the tests are complete, results can be checked by running the
results.csh script in the [suite_name].[testid]::

  cd testsuite.[testid]
  ./results.csh

To report the test results, as is required for Pull Requests to be accepted into 
the master the CICE Consortium code see :ref:`testreporting`.

Multiple suites are supported on the command line as comma separated arguments::

  ./cice.setup --suite base_suite,decomp_suite --mach wolf --env gnu --testid myid

If a user adds ``--set`` to the suite, all tests in that suite will add that option::

  ./cice.setup --suite base_suite,decomp_suite --mach wolf --env gnu --testid myid -s debug

The option settings defined in the suite have precendent over the command line
values if there are conflicts.

The predefined test suites are defined under **configuration/scripts/tests** and 
the files defining the suites
have a suffix of .ts in that directory.  The format for the test suite file 
is relatively simple.  
It is a text file with white space delimited 
columns that define a handful of values in a specific order.  
The first column is the test name, the second the grid, the third the pe count, 
the fourth column is
the ``--set`` options and the fifth column is the ``--diff`` argument. 
The fourth and fifth columns are optional.
Lines that begin with # or are blank are ignored.  For example,
::

   #Test   Grid  PEs  Sets                Diff
    smoke   col  1x1  diag1  
    smoke   col  1x1  diag1,run1year  smoke_col_1x1_diag1
    smoke   col  1x1  debug,run1year  
   restart  col  1x1  debug  
   restart  col  1x1  diag1  
   restart  col  1x1  pondcesm  
   restart  col  1x1  pondlvl  
   restart  col  1x1  pondtopo  

The argument to ``--suite`` defines the test suite (.ts) filename and that argument 
can contain a path.  
**cice.setup** 
will look for the filename in the local directory, in **configuration/scripts/tests/**, 
or in the path defined by the ``--suite`` option.

Because many of the command line options are specified in the input file, ONLY the
following options are valid for suites,

``--suite`` filename
  required, input filename with list of suites

``--mach`` MACHINE
  required

``--env`` ENVIRONMENT1,ENVIRONMENT2
  strongly recommended

``--set`` SET1,SET2
  optional

``--acct`` ACCOUNT
  optional

``--testid`` ID
  required

``--bdir`` DIR
  optional, top level baselines directory and defined by default by ICE_MACHINE_BASELINE in **env.[machine]_[environment]**.

``--bgen`` DIR
  recommended, test output is copied to this directory under [bdir]

``--bcmp`` DIR
  recommended, test output are compared to prior results in this directory under [bdir]

``--report``
  This is only used by ``--suite`` and when set, invokes a script that sends the test results to the results page when all tests are complete.  Please see :ref:`testreporting` for more information.

Please see :ref:`case_options` and :ref:`indtests` for more details about how these options are used.


Test Suite Examples
~~~~~~~~~~~~~~~~~~~~~~~~

 1) **Basic test suite**
     
    Specify suite, mach, env, testid.
    ::

     ./cice.setup --suite base_suite --mach conrad --env cray --testid v01a
     cd base_suite.v01a
     #wait for runs to complete
     ./results.csh

 2) **Basic test suite on multiple environments**

    Specify multiple envs.
    ::

      ./cice.setup --suite base_suite --mach conrad --env cray,pgi,intel,gnu --testid v01a
      cd base_suite.v01a
      #wait for runs to complete
      ./results.csh

    Each env can be run as a separate invokation of `cice.setup` but if that
    approach is taken, it is recommended that different testids be used.

 3) **Basic test suite with generate option defined**

    Add ``--set``
    ::

       ./cice.setup --suite base_suite --mach conrad --env gnu --testid v01b --set diag1
       cd base_suite.v01b
       #wait for runs to complete
      ./results.csh

    If there are conflicts between the ``--set`` options in the suite and on the command line,
    the suite will take precedent.

 4) **Multiple test suites from a single command line**

    Add comma delimited list of suites
    ::

      ./cice.setup --suite base_suite,decomp_suite --mach conrad --env gnu --testid v01c
      cd base_suite.v01c
      #wait for runs to complete
      ./results.csh

     If there are redundant tests in multiple suites, the scripts will understand that and only
     create one test.

 5) **Basic test suite, store baselines in user defined name**

    Add ``--bgen``
    ::

      ./cice.setup --suite base_suite --mach conrad --env cray --testid v01a --bgen cice.v01a
      cd base_suite.v01a
      #wait for runs to complete
      ./results.csh

     This will store the results in the default [bdir] directory under the subdirectory cice.v01a.

 6) **Basic test suite, store baselines in user defined top level directory**

    Add ``--bgen`` and ``--bdir``
    ::

      ./cice.setup --suite base_suite --mach conrad --env cray --testid v01a --bgen cice.v01a --bdir /tmp/user/CICE_BASELINES
      cd base_suite.v01a
      #wait for runs to complete
      ./results.csh

    This will store the results in /tmp/user/CICE_BASELINES/cice.v01a.

 7) **Basic test suite, store baselines in auto-generated directory**

    Add ``--bgen default``
    ::

      ./cice.setup --suite base_suite --mach conrad --env cray --testid v01a --bgen default
      cd base_suite.v01a
      #wait for runs to complete
      ./results.csh

     This will store the results in the default [bdir] directory under a directory name generated by the script
     that includes the hash and date.

 8) **Basic test suite, compare to prior baselines**

    Add ``--bcmp``
    ::

      ./cice.setup --suite base_suite --mach conrad --env cray --testid v02a --bcmp cice.v01a
      cd base_suite.v02a
      #wait for runs to complete
      ./results.csh

    This will compare to results saved in the baseline [bdir] directory under
    the subdirectory cice.v01a. With the ``--bcmp`` option, the results will be tested
    against prior baselines to verify bit-for-bit, which is an important step prior 
    to approval of many (not all, see :ref:`compliance`) Pull Requests to incorporate code into 
    the CICE Consortium master code. You can use other regression options as well.
    (``--bdir`` and ``--bgen``)

 9) **Basic test suite, use of default string in regression testing**

    default is a special argument to ``--bgen`` and ``--bcmp``.  When used, the
    scripts will automate generation of the directories.  In the case of ``--bgen``,
    a unique directory name consisting of the hash and a date will be created.
    In the case of ``--bcmp``, the latest directory in [bdir] will automatically
    be used.  This provides a number of useful features

     - the ``--bgen`` directory will be named after the hash automatically
     - the ``--bcmp`` will always find the most recent set of baselines
     - the ``--bcmp`` reporting will include information about the comparison directory 
       name which will include hash information
     - automation can be invoked easily, especially if ``--bdir`` is used to create separate
       baseline directories as needed.

    Imagine the case where the default settings are used and ``--bdir`` is used to 
    create a unique location.  You could easily carry out regular builds automatically via,
    ::

      set mydate = `date -u "+%Y%m%d"`
      git clone https://github.com/myfork/cice cice.$mydate --recursive
      cd cice.$mydate
      ./cice.setup --suite base_suite --mach conrad --env cray,gnu,intel,pgi --testid $mydate --bcmp default --bgen default --bdir /tmp/work/user/CICE_BASELINES_MASTER

    When this is invoked, a new set of baselines will be generated and compared to the prior
    results each time without having to change the arguments.

 10) **Create and test a custom suite**

    Create your own input text file consisting of 5 columns of data,
     - Test
     - Grid
     - pes
     - sets (optional)
     - diff test (optional)

    such as
    ::

       > cat mysuite
       smoke    col  1x1  diag1,debug
       restart  col  1x1
       restart  col  1x1  diag1,debug    restart_col_1x1
       restart  col  1x1  mynewoption,diag1,debug

    then use that input file, mysuite
    ::

      ./cice.setup --suite mysuite --mach conrad --env cray --testid v01a --bgen default
      cd mysuite.v01a
      #wait for runs to complete
      ./results.csh

    You can use all the standard regression testing options (``--bgen``, ``--bcmp``, 
    ``--bdir``).  Make sure any "diff" testing that goes on is on tests that
    are created earlier in the test list, as early as possible.  Unfortunately,
    there is still no absolute guarantee the tests will be completed in the correct 
    sequence.


.. _testreporting:

Test Reporting
---------------

The CICE testing scripts have the capability to post test results
to the official CICE Consortium Test-Results 
`wiki page <https://github.com/CICE-Consortium/Test-Results/wiki>`_.
You may need write permission on the wiki. If you are interested in using the
wiki, please contact the consortium. Note that in order for code to be 
accepted to the CICE Consortium master through a Pull Request it is necessary
for the developer to provide proof that their code passes relevant tests. 

To post results, once a test suite is complete, run ``results.csh`` and
``report_results.csh`` from the suite directory,
::

  ./cice.setup --suite base_suite --mach conrad --env cray --testid v01a
  cd base_suite.v01a
  #wait for runs to complete
  ./results.csh
  ./report_results.csh

The reporting can also be automated by adding ``--report`` to ``cice.setup``
::

  ./cice.setup --suite base_suite --mach conrad --env cray --testid v01a --report

With ``--report``, the suite will create all the tests, build and submit them,
wait for all runs to be complete, and run the results and report_results scripts.


.. _compliance:

Code Compliance Test (non bit-for-bit validation)
----------------------------------------------------

A core tenet of CICE dycore and CICE innovations is that they must not change 
the physics and biogeochemistry of existing model configurations, notwithstanding 
obsolete model components. Therefore, alterations to existing CICE Consortium code
must only fix demonstrable numerical or scientific inaccuracies or bugs, or be 
necessary to introduce new science into the code.  New physics and biogeochemistry 
introduced into the model must not change model answers when switched off, and in 
that case CICEcore and CICE must reproduce answers bit-for-bit as compared to 
previous simulations with the same namelist configurations. This bit-for-bit 
requirement is common in Earth System Modeling projects, but often cannot be achieved 
in practice because model additions may require changes to existing code.  In this 
circumstance, bit-for-bit reproducibility using one compiler may not be unachievable 
on a different computing platform with a different compiler.  Therefore, tools for 
scientific testing of CICE code changes have been developed to accompany bit-for-bit 
testing. These tools exploit the statistical properties of simulated sea ice thickness 
to confirm or deny the null hypothesis, which is that new additions to the CICE dycore 
and CICE have not significantly altered simulated ice volume using previous model 
configurations.  Here we describe the CICE testing tools, which are applies to output 
from five-year gx-1 simulations that use the standard CICE atmospheric forcing. 
A scientific justification of the testing is provided in
:cite:`Hunke18`. The following sections follow :cite:`Roberts18`.

.. _paired:


Two-Stage Paired Thickness Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first quality check aims to confirm the null hypotheses
:math:`H_0\!:\!\mu_d{=}0` at every model grid point, given the mean
thickness difference :math:`\mu_d` between paired CICE simulations
‘:math:`a`’ and ‘:math:`b`’ that should be identical. :math:`\mu_d` is
approximated as
:math:`\bar{h}_{d}=\tfrac{1}{n}\sum_{i=1}^n (h_{ai}{-}h_{bi})` for
:math:`n` paired samples of ice thickness :math:`h_{ai}` and
:math:`h_{bi}` in each grid cell of the gx-1 mesh. Following
:cite:`Wilks06`, the associated :math:`t`-statistic
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

Following :cite:`Zwiers95`, the effective sample size is
limited to :math:`n_{eff}\in[2,n]`. This definition of :math:`n_{eff}`
assumes ice thickness evolves as an AR(1) process
:cite:`vonstorch99`, which can be justified by analyzing
the spectral density of daily samples of ice thickness from 5-year
records in CICE Consortium member models :cite:`Hunke18`.
The AR(1) approximation is inadmissible for paired velocity samples,
because ice drift possesses periodicity from inertia and tides
:cite:`Hibler06,Lepparanta12,Roberts15`. Conversely,
tests of paired ice concentration samples may be less sensitive to ice
drift than ice thickness. In short, ice thickness is the best variable
for CICE Consortium quality control (QC), and for the test of the mean
in particular.

Care is required in analyzing mean sea ice thickness changes using
(:eq:`t-distribution`) with
:math:`N{=}n_{eff}{-}1` degrees of freedom.
:cite:`Zwiers95` demonstrate that the :math:`t`-test in
(:eq:`t-distribution`) becomes conservative when
:math:`n_{eff} < 30`, meaning that :math:`H_0` may be erroneously
confirmed for highly auto-correlated series. Strong autocorrelation
frequently occurs in modeled sea ice thickness, and :math:`r_1>0.99` is
possible in parts of the gx-1 domain for the five-year QC simulations.
In the event that :math:`H_0` is confirmed but :math:`2\leq n_{eff}<30`,
the :math:`t`-test progresses to the ‘Table Lookup Test’ of
:cite:`Zwiers95`, to check that the first-stage test
using (:eq:`t-distribution`) was not
conservative. The Table Lookup Test chooses critical :math:`t` values
:math:`|t|<t_{crit}({1{-}\alpha/2},N)` at the :math:`\alpha`
significance level based on :math:`r_1`. It uses the conventional
:math:`t={\bar{h}_{d} \sqrt{n}}/{\sigma_d}` statistic with degrees of
freedom :math:`N{=}n{-}1`, but with :math:`t_{crit}` values generated
using the Monte Carlo technique described in
:cite:`Zwiers95`, and summarized in :ref:`Table-Lookup` for 5-year QC
simulations (:math:`N=1824`) at the two-sided 80% confidence interval
(:math:`\alpha=0.2`). We choose this interval to limit Type II errors,
whereby a QC test erroneously confirms :math:`H_0`.

Table :ref:`Table-Lookup` shows the summary of two-sided :math:`t_{crit}` values for the Table
Lookup Test of :cite:`Zwiers95` at the 80% confidence
interval generated for :math:`N=1824` degrees of freedom and lag-1
autocorrelation :math:`r_1`.

.. _Table-Lookup:

.. csv-table:: Two-sided :math:`t_{crit}` values
   :widths: 10, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5

   :math:`r_1`,-0.05,0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.97,0.99
   :math:`t_{crit}`,1.32,1.32,1.54,2.02,2.29,2.46,3.17,3.99,5.59,8.44,10.85,20.44


.. _quadratic:


Quadratic Skill Compliance Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the two-stage test of mean sea ice thickness, we also
check that paired simulations are highly correlated and have similar
variance using a skill metric adapted from
:cite:`Taylor01`. A general skill score applicable to
Taylor diagrams takes the form

.. math::
   S_m=\frac{4(1+R)^m}{({\hat{\sigma}_{f}+1/{\hat{\sigma}_{f}}})^2 (1+R_0)^m}
   :label: taylor-skill

where :math:`m=1` for variance-weighted skill, and :math:`m=4` for
correlation-weighted performance, as given in equations (4) and (5) of
:cite:`Taylor01`, respectively. We choose :math:`m=2` to
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
provided in :cite:`Hunke18`.


Code Compliance Testing Procedure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CICE code compliance test is performed by running a python script 
(**configurations/scripts/tests/QC/cice.t-test.py**).
In order to run the script, the following requirements must be met:

* Python v2.7 or later
* netCDF Python package
* numpy Python package
* matplotlib Python package (optional)
* basemap Python package (optional)

In order to generate the files necessary for the compliance test, test cases should be
created with the ``qc`` option (i.e., ``--set qc``) when running cice.setup.  This 
option results in daily, non-averaged history files being written for a 5 year simulation.

To install the necessary Python packages, the ``pip`` Python utility can be used.

.. code-block:: bash

  pip install --user netCDF4
  pip install --user numpy
  pip install --user matplotlib

To run the compliance test, setup a baseline run with the original baseline model and then 
a perturbation run based on recent model changes.  Use ``--sets qc`` in both runs in addition
to other settings needed.  Then use the QC script to compare history output,

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

The cice.t-test.py requires memory to store multiple two-dimensional fields spanning 
1825 unique timesteps, a total of several GB.  An appropriate resource is needed to 
run the script.  If the script runs out of memory on an interactive resource, try
logging into a batch resource or finding a large memory node.


End-To-End Testing Procedure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example of a step-by-step procedure for testing a code change that might result in non bit-for-bit results.   First, run a regression test,

.. code-block:: bash

  # Run a full regression test to verify bit-for-bit

  # Create a baseline dataset (only necessary if no baseline exists on the system)
  # git clone the baseline code

  ./cice.setup -m onyx -e intel --suite base_suite --testid base0 -bgen cice.my.baseline

  # Run the test suite with the new code
  # git clone the new code

  ./cice.setup -m onyx -e intel --suite base_suite --testid test0 --bcmp cice.my.baseline

  # Check the results

  cd testsuite.test0
  ./results.csh

..

If the regression comparisons fail, then you may want to run the QC test,

.. code-block:: bash

  # Run the QC test

  # Create a QC baseline
  # From the baseline sandbox

  ./cice.setup -m onyx -e intel --test smoke -g gx1 -p 44x1 --testid qc_base -s qc,medium
  cd onyx_intel_smoke_gx1_44x1_medium_qc.qc_base
  ./cice.build
  ./cice.submit

  # Create the t-test testing data
  # From the update sandbox

  ./cice.setup -m onyx -e intel --test smoke -g gx1 -p 44x1 -testid qc_test -s qc,medium
  cd onyx_intel_smoke_gx1_44x1_medium_qc.qc_test
  ./cice.build
  ./cice.submit

  # Wait for runs to finish
  # Perform the QC test

  cp configuration/scripts/tests/QC/cice.t-test.py
  ./cice.t-test.py /p/work/turner/CICE_RUNS/onyx_intel_smoke_gx1_44x1_medium_qc.qc_base \
                   /p/work/turner/CICE_RUNS/onyx_intel_smoke_gx1_44x1_medium_qc.qc_test

  # Example output:
  INFO:__main__:Number of files: 1825
  INFO:__main__:Two-Stage Test Passed
  INFO:__main__:Quadratic Skill Test Passed for Northern Hemisphere
  INFO:__main__:Quadratic Skill Test Passed for Southern Hemisphere
  INFO:__main__:
  INFO:__main__:Quality Control Test PASSED


.. _testplotting:

Test Plotting
----------------

The CICE scripts include a script (``timeseries.csh``) that will generate a timeseries 
figure from the diagnostic output file.  
When running a test suite, the ``timeseries.csh`` script is automatically copied to the suite directory.  
If the ``timeseries.csh`` script is to be used on a test / case that is not a part of a test suite, 
users will need to run the ``timeseries.csh`` script from the tests directory 
(``./configuration/scripts/tests/timeseries.csh``), or copy it to a local directory and run it 
locally (``cp configuration/scripts/tests/timeseries.csh .`` followed by 
``./timeseries.csh /path/to/ice_diag.full_ITD``. The plotting script can be run
on any of the output files - icefree, slab, full_ITD, land).  To generate the figure, 
run the ``timeseries.csh`` script and pass the full path to the ice_diag file as an argument.  

For example:

Run the test suite. ::

$ ./cice.setup -m conrad -e intel --suite base_suite --testid t00

Wait for suite to finish then go to the directory. ::

$ cd base_suite.t00

Run the timeseries script on the desired case. ::

$ ./timeseries.csh /p/work1/turner/CICE_RUNS/conrad_intel_smoke_col_1x1_diag1_run1year.t00/ice_diag.full_ITD
    
The output figures are placed in the directory where the ice_diag file is located.

This plotting script can be used to plot the following variables:

  - area fraction
  - average ice thickness (m)
  - average snow depth (m)
  - air temperature (C)
  - shortwave radiation (:math:`W/m^2`)
  - longwave radiation (:math:`W/m^2`)
  - snowfall
  - average salinity (ppt)
  - surface temperature (C)
  - outward longwave flux (:math:`W/m^2`)
  - sensible heat flux (:math:`W/m^2`)
  - latent heat flux (:math:`W/m^2`)
  - top melt (m)
  - bottom melt (m)
  - lateral melt (m)
  - new ice (m)
  - congelation (m)
  - snow-ice (m)
  - initial energy change (:math:`W/m^2`)

