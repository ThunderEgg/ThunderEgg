# thunderegg
![alt text](https://github.com/GEM3D/pressurePoissonSolver/blob/master/icon.png)

ThunderEgg is an object-oriented C++ library designed for flexibility and to allow users to implement parallel multigrid preconditioners for various problems on octree and quadtree adaptive meshes.

# Members of the team :

* Scott Aiton
* Donna Calhoun
* Grady Wright

# Required Software
* MPI
* CMake
* BLAS and LAPACK

## Optional Software
* FFTW - [FFTW PatchSolver](https://thunderegg.dev/docs/develop/doc/html/classThunderEgg_1_1Poisson_1_1FFTWPatchSolver.html)
* PETSc - ThunderEgg provides [a set of interfaces](https://thunderegg.dev/docs/develop/doc/html/namespaceThunderEgg_1_1PETSc.html) to use PETSc Krylov Solvers and PETSc matrices.
* p4est - for compatibility with the p4est quadtree library

# Compiling
Create a seperate source directory and run cmake in the build directory:
```
$ cd build_dir
$ cmake /path/to/source
```
Compilers can be specified in the following way
```
$ cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER /path/to/source
```
Paths to libraries can be specified with the following cmake options
```
-DFFTW_DIR=/path/to/library
-DPETSC_DIR=/path/to/library
-Dp4est_DIR=/path/to/library

```
Then compile with make:
```
make
```

# Examples:

Example applications are in the `apps` directory.

TODO cleanup/simplify example applications