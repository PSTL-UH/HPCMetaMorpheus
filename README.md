HPCMorpheus

This is a C++ version of the MetaMorpheus library, targetting HPC environments. 

The original C# version was developed by the L.M. Smith Group https://github.com/smith-chem-wisc/
This is currently work in progress. The code has been translated to C++ by the PSTL group at the University of Houston, 
and will be parallelized in the future. The software is distributed under a MIT license.

The repository contains 2 branches, the main (master) sequential branch with the C++ code, and the multi-threaded OpenMP based source (topic/OpenMP-version). 
An MPI version for clusters is currently being conceptually developed.

As of right now (September 2020), only XLSeaarchTask has received significant testing and considered reasonably stable, although bugs and issues most likely still exist. 
SearchTask and CalibrationTasks are work in progress should be available in the near future.

