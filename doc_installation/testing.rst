

Testing the installation
========================

It is very important that you verify that your Dalton installation correctly
reproduces the reference test set before running any production calculations.

The test set driver is CTest which can be invoked with "make test" after building
the code.


Environment variables for testing
---------------------------------

Before testing with "make test" you should export the
following environment variables::

  $ export DALTON_TMPDIR=/scratch        # scratch space for Dalton (adapt the path of course)
  $ export CTEST_MAKE_NUM_PROCS=16       # in this case the code will be compiled with 16 processes (make -j16)
  $ export DALTON_NUM_MPI_PROCS=4        # in this case 4 processes, only relevant if you compile with MPI

Note that if you set the DALTON_NUM_MPI_PROCS to something different from 1,
the dalton script will assume you have compiled using MPI and run the mpirun
command!


Running the test set
--------------------

You can run the whole test set either using::

  $ make test

or directly through CTest::

  $ ctest

Both are equivalent ("make test" runs CTest) but running
CTest directly makes it easier to run sequential tests on several
cores::

  $ ctest -j4

You can select the subset of tests by matching test names to a regular expression::

  $ ctest -R dft

Alternatively you can select the tests with a label matching a regular expression::

  $ ctest -L rsp

The following command will give you all available labels::

  $ ctest --print-labels


Running only DALTON or only LSDALTON tests
------------------------------------------

Only DALTON tests::

  $ ctest -L dalton

Only LSDALTON tests::

  $ ctest -L linsca

All tests::

  $ ctest
