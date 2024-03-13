if(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)

    target_sources(unit_tests_mpi1 PRIVATE DFTPatchSolver_MPI1.cpp)

endif(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)

if(TARGET PETSc::PETSc AND TARGET FFTW::FFTW)

    target_sources(unit_tests_mpi1 PRIVATE FastSchurMatrixAssemble2D_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE FastSchurMatrixAssemble3D_MPI1.cpp)

endif(TARGET PETSc::PETSc AND TARGET FFTW::FFTW)

if(TARGET FFTW::FFTW)

    target_sources(unit_tests_mpi1 PRIVATE FFTWPatchSolver_MPI1.cpp)

endif(TARGET FFTW::FFTW)

if(TARGET PETSc::PETSc)

    target_sources(unit_tests_mpi1 PRIVATE MatrixHelper_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE MatrixHelper2d_MPI1.cpp)

endif(TARGET PETSc::PETSc)

target_sources(unit_tests_mpi1 PRIVATE StarPatchOperator_MPI1.cpp)