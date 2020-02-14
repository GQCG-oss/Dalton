Nightly runs on full testset: https://testboard.org/cdash/index.php?project=Dalton

## Dalton links

- [Home page](http://daltonprogram.org/)
- [Forum](http://forum.daltonprogram.org/)
- [Article](http://onlinelibrary.wiley.com/doi/10.1002/wcms.1172/abstract)


## Quick start

Note that it is currently not practical to download the source using the
download button on GitLab, because it will not include the submodules that are
required to build Dalton. Instead you should clone the repository as described
below.

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
where you replace *Dalton2018.0* by the release version that you are
interested in. The list of past releases available in this repository can be
found here: https://gitlab.com/dalton/dalton/-/tags.

In case you did not include the `--recursive` argument when you cloned the
repository, it is necessary to run the following two commands:
```
$ git submodule init
$ git submodule update --recursive
```

To build the code, perform the following steps:
```
$ ./setup [--help]
$ cd build
$ make [-j4]
```

Once the build is complete, you can run the test set as:
```
$ ctest [-j4]
```

To switch branch (or release tag), run the following two commands from the *dalton* directory:
```
$ git checkout feature-branch
$ git submodule update
```
This can also be achieved in one step when you clone the repository:
```
$ git clone --recursive -b feature-branch https://gitlab.com/dalton/dalton.git
```
