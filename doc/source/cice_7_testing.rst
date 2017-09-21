..

Testing CICE
================

Version 6, August 2017
This documents how to use the testing features developed for the 
CICE Consortium CICE sea ice model.

.. _basic:

Individual tests and test suites
--------------------------------

The CICE scripts support both setup of individual tests as well as test suites.  Individual
tests are run from the command line like

  > create.case -t smoke -m wolf -g gx3 -p 8x2 -s diag1,run5day -testid myid

where -m designates a specific machine.  Test suites are multiple tests that are specified in 
an input file and are started on the command line like

  > create.case -ts base_suite -m wolf -testid myid

create.case with -t or -ts require a testid to uniquely name test directories.  The format
of the case directory name for a test will always be 
${machine}_${test}_${grid}_${pes}_${soptions}.${testid}

To build and run a test, the process is the same as a case,
  cd into the test directory,
  
  run cice.build
  
  run cice.submit

The test results will be generated in a local file called "test_output".

When running a test suite, the create.case command line automatically generates all the tests
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

Please run './create.case -h' for additional details.

.. _additional:

Additional testing options
--------------------------

There are several additional options on the create.case command line for testing that
provide the ability to regression test and compare tests to each other.

  -bd defines a baseline directory where tests can be stored for regression testing
  
  -bg defines a version name that where the current tests can be saved for regression testing
  
  -bc defines a version name that the current tests should be compared to for regression testing
  
  -td provides a way to compare tests with each other

To use -bg,
  > create.case -ts base_suite -m wolf -testid v1 -bg version1 -bd $SCRATCH/CICE_BASELINES
will copy all the results from the test suite to $SCRATCH/CICE_BASELINES/version1.

To use -bc,
  > create.case -ts base_suite -m wolf -testid v2 -bc version1 -bd $SCRATCH/CICE_BASELINES
will compare all the results from this test suite to results saved before in $SCRATCH/CICE_BASELINES/version1.

-bc and -bg can be combined,
  >create.case -ts base_suite -m wolf -testid v2 -bg version2 -bc version1 -bd $SCRATCH/CICE_BASELINES
will save the current results to $SCRATCH/CICE_BASELINES/version2 and compare the current results to
results save before in $SCRATCH/CICE_BASELINES/version1.

-bg, -bc, and -bd are used for regression testing.  There is a default -bd on each machine.

-td allows a user to compare one test result to another.  For instance,
  > create.case -t smoke -m wolf -g gx3 -p 8x2 -s run5day -testid t01
  > create.case -t smoke -m wolf -g gx3 -p 4x2 -s run5day -testid t01 -td smoke_gx3_8x2_run5day
An additional check will be done for the second test (because of the -td argument), and it will compare
the output from the first test "smoke_gx3_8x2_run5day" to the output from it's test "smoke_gx3_4x2_run5day"
and generate a result for that.  It's important that the first test complete before the second test is 
done.  Also, the -td option works only if the testid and the machine are the same for the baseline
run and the current run.

.. _format:

Test suite format
-----------------

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
The argument to -ts defines which filename to choose and that argument can contain a path.  create.case 
will also look for the filename in configuration/scripts/tests where some preset test suites are defined.

Example Tests (Quickstart)
--------------------------

To generate a baseline dataset for a test case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

./create.case -t smoke -m wolf -bg cicev6.0.0 -testid t00

cd wolf_smoke_gx3_4x1.t00

./cice.build

./cice.submit

# After job finishes, check output

cat test_output


To run a test case and compare to a baseline dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

./create.case -t smoke -m wolf -bc cicev6.0.0 -testid t01

cd wolf_smoke_gx3_4x1.t01

./cice.build

./cice.submit

# After job finishes, check output

cat test_output


To run a test suite to generate baseline data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

./create.case -m wolf -ts base_suite -testid t02 -bg cicev6.0.0bs

cd base_suite.t02

# Once all jobs finish, concatenate all output

./results.csh  # All tests results will be stored in results.log

# To plot a timeseries of "total ice extent", "total ice area", and "total ice volume"

./timeseries.csh <directory>

ls *.png


To run a test suite to compare to baseline data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

./create.case -m wolf -ts base_suite -testid t03 -bc cicev6.0.0bs

cd base_suite.t03

# Once all jobs finish, concatenate all output

./results.csh  # All tests results will be stored in results.log

# To plot a timeseries of "total ice extent", "total ice area", and "total ice volume"

./timeseries.csh <directory>

ls *.png


To compare to another test
~~~~~~~~~~~~~~~~~~~~~~~~~~
`First:`

./create.case -m wolf -t smoke -testid t01 -p 8x2

cd wolf_smoke_gx3_8x2.t01

./cice.build

./cice.submit

# After job finishes, check output

cat test_output

`Then, do the comparison:` 

./create.case -m wolf -t smoke -testid t01 -td smoke_gx3_8x2 -s thread -p 4x1

cd wolf_smoke_gx3_4x1_thread.t01

./cice.build

./cice.submit

# After job finishes, check output

cat test_output


Additional Details
------------------
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
    end of the case directory.  For example, "./create.case -m wolf -t smoke -testid t12 -p 4x1"
    creates the directory wolf_smoke_gx3_4x1.t12.  This flag is REQUIRED if using -t or -ts.
