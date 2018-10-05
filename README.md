[![pipeline status](https://gitlab.com/dalton/dalton/badges/master/pipeline.svg)](https://gitlab.com/dalton/dalton/pipelines) [![coverage report](https://gitlab.com/dalton/dalton/badges/master/coverage.svg?job=nightly-coverage)](https://gitlab.com/dalton/dalton/pipelines)

Nightly runs on full testset: https://testboard.org/cdash/index.php?project=Dalton

# Quick start

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

# How to contribute

See Dalton Developer’s Guide: http://dalton-devguide.readthedocs.io

# Dalton links

- [Home page](http://daltonprogram.org/)
- [Forum](http://forum.daltonprogram.org/)
- [Article](http://onlinelibrary.wiley.com/doi/10.1002/wcms.1172/abstract)
