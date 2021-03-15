list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/StarPatchOperator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/StarPatchOperator.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/DFTPatchSolver.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/DFTPatchSolver.cpp)

if(FFTW_FOUND)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/FFTWPatchSolver.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/FFTWPatchSolver.cpp)

endif(FFTW_FOUND)

if(PETSC_FOUND)
  if(FFTW_FOUND)
    list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/FastSchurMatrixAssemble2D.h)
    list(APPEND ThunderEgg_SRCS
                ThunderEgg/Poisson/FastSchurMatrixAssemble2D.cpp)

    list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/FastSchurMatrixAssemble3D.h)
    list(APPEND ThunderEgg_SRCS
                ThunderEgg/Poisson/FastSchurMatrixAssemble3D.cpp)
  endif(FFTW_FOUND)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/MatrixHelper.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/MatrixHelper.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/MatrixHelper2d.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/MatrixHelper2d.cpp)

endif(PETSC_FOUND)
