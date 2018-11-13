[![pipeline status](https://gitlab.com/dalton/dalton/badges/master/pipeline.svg)](https://gitlab.com/dalton/dalton/pipelines)
[![coverage report](https://gitlab.com/dalton/dalton/badges/master/coverage.svg?job=nightly-coverage)](https://gitlab.com/dalton/dalton/pipelines)
[![codecov](https://codecov.io/gl/dalton/dalton/branch/master/graph/badge.svg)](https://codecov.io/gl/dalton/dalton)

Nightly runs on full testset: https://testboard.org/cdash/index.php?project=Dalton

## Quick start

Clone the repository:
```
$ git clone --recursive https://gitlab.com/dalton/dalton.git
```
Build the code:
```
$ cd dalton
$ ./setup [--help]
$ cd build
$ make [-j4]
```

Run the test set:
```
$ ctest [-j4]
```

To switch branch to, e.g., the latest release run the following two commands from the *dalton* directory:
```
$ git checkout release/2016
$ git submodule update
```
This can also be achieved in one step when you clone the repository:
```
$ git clone --recursive -b release/2016 https://gitlab.com/dalton/dalton.git
```
In case you did not include the `--recursive` argument when you cloned the repository, it is necessary to run the following two commands from the *dalton* directory before entering the *build* directory and running `make`:
```
$ git submodule init
$ git submodule update --recursive
```

## How to contribute

See Dalton Developer’s Guide: http://dalton-devguide.readthedocs.io

## Dalton links

- [Home page](http://daltonprogram.org/)
- [Forum](http://forum.daltonprogram.org/)
- [Article](http://onlinelibrary.wiley.com/doi/10.1002/wcms.1172/abstract)
