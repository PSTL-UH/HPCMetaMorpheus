Compiling and installing the HPCMetaMorpheus library follows standard Linux software practices, requiring a configure and a make step.

The configure step will check for all required dependencies, which are as of right now the following libraries:

1. A C++ compiler support the C++2017 standard or newer (e.g. gcc 8.1 or later) 
2. boost math library
3. xml2 library
4. xsdcxx and xerces

An example configuration line could look as follows:

```
./configure --disable-debug
make -j 8
```

The software has been tested with gcc 8.2, gcc9.3, gcc 10.2, and intel icpc 19.0.5 compilers, on multiple Linux systems (Redhat 7.6, OpenSuSe 15.2, Arch Linux, Ubuntu). Note: since this is still work in progres, the library is compiled by default with debug symbols and without optimization. For any form of performance evaluation, the user has to disable these setting by using the ```--disable-debug``` option during the configure step. Configuration options and default values might change in the future.

The HPCMzlib library is incorporated as a git submodule in HPCMetaMorpheus. There are a few items to consider when using git submodules, the most important ones being that is not sufficient to simply `git clone <GIT_URL_REPO>` a repository.  Instead, you must add the `--recursive` option to your clone command, e.g.:

```shell
$ git clone --recursive git@github.com:PSTL-UH/HPCMetaMorpheus.git
```

If you already cloned the HPCMetaMorpheus repository and didn't use `--recursive`, you can also download submodules retroactively, e.g.:

```shell
$ git clone  git@github.com:PSTL-UH/HPCMetaMorpheus.git
$ git submodule update --init --recursive
```

Note that `git pull ...` is also no longer sufficient to update your source tree.  This will still update everything inside the HPCMetaMorpheus repository, but it will not (by default) update any changes to submodules. The easiest solution to fetch all changes that have been applied to submodules as well is to use `--recurse-submodules` option, e.g.:

```shell
$ git pull --recurse-submodules
```

To checkout and compile the MPI version, please checkout the ```topic/MPI``` branch and provide an MPI capable compile environment. For example:
```
$ git checkout topic/MPI
$ ./configure CXX=mpicxx --disable-debug
```

Similarly, for the OpenMP version checkout the ```topic/OpenMP``` branch, and provide an OpenMP capable compiler.
