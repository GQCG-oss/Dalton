

Environment variables
---------------------

Before starting testing ("make test" or nightly), you should export the
following environment variables::

  export DALTON_TMPDIR=/scratch     # scratch space for Dalton (adapt the path of course)
  export DALTON_NR_MPI_PROCS=4      # in this case 4 cores, only relevant if you compile with MPI


Testing
-------

You can run the test set either using::

  $ make test

or directly through ctest::

  $ ctest

Both are equivalent ("make test" runs ctest) but running
ctest directly makes it easier to run sequential tests on several
cores::

  $ ctest -j4

You can run a subset of tests matching test names to a regular expression::

  $ ctest -R dft*

Alternatively you can run a subtest with a specific label::

  $ ctest -L linsca

The following command will give you all available labels::

  $ ctest --print-labels


Nightly testing
---------------

radovan will write this ...   ctest -D Nightly/Experimental -j4
