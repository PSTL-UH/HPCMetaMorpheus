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

The software has been tested with gcc 8.2, gcc 10.2, and intel icpc 19.0.5, on multiple Linux systems (Redhat 7.6, OpenSuSe 15.2, Arch Linux)
Note: since this is still work in development, the library is by default compiled with debug symbols and without optimization. For any form of performance evaluation,
the user has to disable the setting by using the ```--disable-debug``` option during the configure step. Configuration values and options might change in the future.

The HPCMzlib library is incorporated as a git submodule in HPCMetaMorpheus. There are a few differences in using submodules, the most important ones being that is not sufficient to simply `git clone <GIT_URL_REPO>`. Instead, you must add `--recursive` into your clone command, e.g.:

```shell
$ git clone --recursive git@github.com:PSTL-UH/HPCMetaMorpheus.git
```

If you already cloned the HPCMetaMorpheus repository and didn't use `--recursive`, you can initialze / download all submodules:

```shell
$ git submodule update --init --recursive
```

Note that `git pull ...` is also no longer sufficient to update your entire source tree.  This will still update everything inside the HPCMetaMorpheus repository, but it will not (by default) update any changes to submodules. The easiest solution is to use `--recurse-submodules`, e.g.:

```shell
$ git pull --recurse-submodules
```
