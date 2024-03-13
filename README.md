# thunderegg

![alt text](https://github.com/GEM3D/pressurePoissonSolver/blob/master/icon.png)

ThunderEgg is an object-oriented C++ library designed for flexibility and to allow users to implement parallel multigrid preconditioners for various problems on octree and quadtree adaptive meshes.

## Members of the team

* Scott Aiton
* Donna Calhoun
* Grady Wright

## Required Software

* MPI
* CMake

## Optional Software
* FFTW - [FFTW PatchSolver](https://thunderegg.dev/ThunderEgg/docs/develop-wip/classThunderEgg_1_1Poisson_1_1FFTWPatchSolver.html)
* BLAS and LAPACK - [DFT PatchSolver](https://thunderegg.dev/ThunderEgg/docs/develop-wip/classThunderEgg_1_1Poisson_1_1DFTPatchSolver.html)
* PETSc - ThunderEgg provides [a set of interfaces](https://thunderegg.dev/ThunderEgg/docs/develop-wip/namespaceThunderEgg_1_1PETSc.html) to use PETSc Krylov Solvers and PETSc matrices.
* p4est - for compatibility with the p4est quadtree library

## Compiling

Configure CMake project:

```sh
cmake -S /path/to/source -B build
```

Compilers can be specified in the following way

```sh
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -S /path/to/source -B build
```

Then compile:

```sh
cmake --build build
```

Some helpful CMake variables for configuration:

Variable           | Default Value |   Description
-------------------|:---------:|---------------------------------------
PETSC_DIR          |         |    The PETSc directory
PETSC_ARCH         |         |    The PETSc arch
petsc              |   ON    |    allow the use the use of PETSc
petsc_required     |   OFF   |    fail if PETSc is not found
FFTW_ROOT          |         |    The fftw directory
fftw               |   ON    |    allow the use the use of FFTW
fftw_required      |   OFF   |    fail if FFTW is not found
P4EST_ROOT         |         |    The p4est directory
p4est              |   ON    |    allow the use the use of p4est
p4est_required     |   OFF   |    fail if p4est is not found
lapack             |   ON    |    allow the use the use of lapack/blas
lapack_required    |   OFF   |    fail if lapack/blas is not found
