

Installation instructions for system administrators
===================================================

Please read the other installation sections for details
but the installation procedure is basically this::

  $ ./setup [--flags] --prefix=/full/install/path/
  $ cd build
  $ make [-jN]
  $ ctest [-jN]
  $ make install

This will install binaries, run scripts, the basis set library,
as well as tools into the install path.

Advise users to always set a suitable scratch directory::

  $ export DALTON_TMPDIR=/full/path/scratch
