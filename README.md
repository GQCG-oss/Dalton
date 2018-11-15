[![pipeline status](https://gitlab.com/dalton/dalton/badges/master/pipeline.svg)](https://gitlab.com/dalton/dalton/pipelines)
[![coverage report](https://gitlab.com/dalton/dalton/badges/master/coverage.svg?job=nightly-coverage)](https://gitlab.com/dalton/dalton/pipelines)
[![codecov](https://codecov.io/gl/dalton/dalton/branch/master/graph/badge.svg)](https://codecov.io/gl/dalton/dalton)

Nightly runs on full testset: https://testboard.org/cdash/index.php?project=Dalton

## Quick start

Clone the repository:
```
$ git clone --recursive https://gitlab.com/dalton/dalton.git
```

This will fetch the entire repository in a directory called *dalton*. By default
it checks out the master branch which is the main development branch. To
checkout a specific release version, run the following commands from inside the
*dalton* directory:
```
$ git checkout Dalton2018.0
$ git submodule update
```
replacing *Dalton2018.0* by the release version that you are interesed in. Note
that it is currently not possible to fetch all past releases.

To build the code:
```
$ ./setup [--help]
$ cd build
$ make [-j4]
```

Run the test set:
```
$ ctest [-j4]
```

To switch branch, run the following two commands from the *dalton* directory:
```
$ git checkout feature-branch
$ git submodule update
```
This can also be achieved in one step when you clone the repository:
```
$ git clone --recursive -b feature-branch https://gitlab.com/dalton/dalton.git
```
In case you did not include the `--recursive` argument when you cloned the repository, it is necessary to run the following two commands from the *dalton* directory before entering the *build* directory and running `make`:
```
$ git submodule init
$ git submodule update --recursive
```

Note that it is currently not practical to download the source using the
download button on GitLab, because it will not include the submodules that are
required to build Dalton. Instead you should clone the repository as described
above.


## How to contribute

See Dalton Developer’s Guide: http://dalton-devguide.readthedocs.io

## Dalton links

- [Home page](http://daltonprogram.org/)
- [Forum](http://forum.daltonprogram.org/)
- [Article](http://onlinelibrary.wiley.com/doi/10.1002/wcms.1172/abstract)
