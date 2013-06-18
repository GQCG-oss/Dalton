

Environment variables
---------------------

Before starting testing ("make test" or nightly), you should export the
following environment variables::

  $ export DALTON_TMPDIR=/scratch        # scratch space for Dalton (adapt the path of course)
  $ export DALTON_NR_MPI_PROCS=4         # in this case 4 cores, only relevant if you compile with MPI

Note that if you set the DALTON_NR_MPI_PROCS to something different from 1 it will assume you have compiled using MPI and run the mpirun command!

Testing
-------

You can run the test set either using::

  $ make test

or directly through CTest::

  $ ctest

Both are equivalent ("make test" runs CTest) but running
CTest directly makes it easier to run sequential tests on several
cores::

  $ ctest -j4

You can run a subset of tests matching test names to a regular expression::

  $ ctest -R dft*

Alternatively you can run a subtest with a specific label::

  $ ctest -L linsca

The following command will give you all available labels::

  $ ctest --print-labels


Running only DALTON or only LSDALTON tests
------------------------------------------

Only DALTON tests::

  $ ctest -L dalton

Only LSDALTON tests::

  $ ctest -L lsdalton

All tests::

  $ ctest


Reporting to the testing dashboard
----------------------------------

We have two testing dashboards, http://repo.ctcc.no/CDash/?project=DALTON, and
http://repo.ctcc.no/CDash/?project=LSDALTON.

By default CTest will report to http://repo.ctcc.no/CDash/?project=DALTON. You
can change this by setting ``CTEST_PROJECT_NAME`` to "LSDALTON" (or some other dashboard)::

  $ export CTEST_PROJECT_NAME=LSDALTON

By default the build name that appears on the dashboard is set to::

  "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${BLAS_TYPE}-${CMAKE_BUILD_TYPE}"

If you don't like it you can either change the default, or set the build name
explicitly::

  $ ./setup -D BUILDNAME='a-better-build-name'

Then run CTest with -D Nightly or Experimental::

  $ ctest -D Nightly      [-jN] [-L ...] [-R ...]
  $ ctest -D Experimental [-jN] [-L ...] [-R ...]

If you want to test your current code, take Experimental. If you want to set
up a cron script to run tests every night, take Nightly.

Currently the above recipes for Nightly/Experimental will compile all targets.
If you want to compile only LSDALTON within Nightly/Experimental, you can do this::

  $ make lsdalton.x
  $ ctest -L lsdalton -D NightlyTest
  $ ctest -L lsdalton -D NightlySubmit

This does not report configure/build problems. We will streamline this soon.
