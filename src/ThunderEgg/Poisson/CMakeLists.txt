list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/StarPatchOperator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/StarPatchOperator.cpp)

if(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)
  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/DFTPatchSolver.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/DFTPatchSolver.cpp)
endif(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)

if(TARGET FFTW::FFTW)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/FFTWPatchSolver.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/FFTWPatchSolver.cpp)

endif(TARGET FFTW::FFTW)

if(PETSC_FOUND)
  if(TARGET FFTW::FFTW)
    list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/FastSchurMatrixAssemble2D.h)
    list(APPEND ThunderEgg_SRCS
                ThunderEgg/Poisson/FastSchurMatrixAssemble2D.cpp)

    list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/FastSchurMatrixAssemble3D.h)
    list(APPEND ThunderEgg_SRCS
                ThunderEgg/Poisson/FastSchurMatrixAssemble3D.cpp)
  endif(TARGET FFTW::FFTW)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/MatrixHelper.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/MatrixHelper.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/MatrixHelper2d.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/MatrixHelper2d.cpp)

endif(PETSC_FOUND)
