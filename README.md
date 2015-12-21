# Quick start

Clone the repository:
```
$ git clone --recursive git@gitlab.com:dalton/dalton.git
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

Dalton Developer’s Guide: http://dalton-devguide.readthedocs.org

## Release tarball creation:

First make sure that you have no local modifications
and no generated files present. For this either clone fresh
or do a careful `rm -rf *; git checkout .; git submodule update --init --recursive`.

Then switch to the release branch.

Finally create the tarball:
```
$ mkdir build
$ cd build
$ cmake ..
$ make release
```