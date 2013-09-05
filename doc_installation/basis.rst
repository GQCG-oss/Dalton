

Basis set directory
===================

The basis set directory is copied to the build directory
(and possibly install directory).
The ``dalton`` and ``lsdalton`` scripts will automatically
find them. You can define or append custom basis set directories
by exporting BASDIR::

  export BASDIR='/somepath:/otherpath'
