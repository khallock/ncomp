NComp (libncomp)
========

"NComp" is a C and Fortran library containing NCL's compiled computational algorithms. NComp provides a C library `libncomp` that can be linked against from Python or any other compiled language with the ability to include C headers.


Build instructions
==================

NComp requires the following dependencies to be installed:

* Any C compiler (GCC and Clang have been tested)
* gfortran
* GNU Autotools/autoconf/automake/libtool (optional)
* CMake (optional, work in progress)

NComp can be built by running the following commands from the root directory of this repository:
```
autoreconf --install
./configure --prefix=$PREFIX
make install
```

where $PREFIX is the path that `ncomp` should be installed.
