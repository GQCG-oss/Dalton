

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
  $ ctest -L linsca -D NightlyTest
  $ ctest -L linsca -D NightlySubmit

This does not report configure/build problems. We will streamline this soon.
