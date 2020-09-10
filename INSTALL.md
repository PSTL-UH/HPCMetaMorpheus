Compiling and installing the HPCMetaMorpheus library follows standard Linux software practices, requiring a configure and a make step.

The configure step will check for all required dependencies, which are as of right now the following libraries:

1. A C++ compiler support the C++2017 standard or newer (e.g. gcc 8.1 or later) 
2. HPCmzlib (both source code and compiled libraries required at the moment).
3. boost math library
4. xml2 library
5. xsdcxx and xerces

An example configuration line could look as follows:

```
./configure --with-mzlib-dir=$(HOME)/HPCmzlib/ --disable-debug
make -j 8
```

The software has been tested with gcc 8.2, gcc 10.2, and intel icpc 19.0.5, on multiple Linux systems (Redhat 7.6, OpenSuSe 15.2, Arch Linux)
Note: since this is still work in development, the library is by default compiled with debug symbols and without optimization. For any form of performance evaluation,
the user has to disable the setting by using the ```--disable-debug``` option during the configure step. Configuration values and options might change in the future.
