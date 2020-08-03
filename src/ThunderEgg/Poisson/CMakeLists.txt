list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/StarPatchOperator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/StarPatchOperator.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/DFTPatchSolver.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/DFTPatchSolver.cpp)

if(FFTW_FOUND)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/FFTWPatchSolver.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/FFTWPatchSolver.cpp)

endif(FFTW_FOUND)

if(PETSC_FOUND)
  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/MatrixHelper.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/MatrixHelper.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/MatrixHelper2d.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/MatrixHelper2d.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/StencilHelper.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/StencilHelper2d.h)
endif(PETSC_FOUND)

list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/Schur/DftPatchSolver.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/Schur/DftPatchSolver.cpp)

if(FFTW_FOUND)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/Schur/FftwPatchSolver.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/Schur/FftwPatchSolver.cpp)

endif(FFTW_FOUND)

list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/Schur/SevenPtPatchOperator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/Schur/SevenPtPatchOperator.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Poisson/Schur/StarPatchOperator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Poisson/Schur/StarPatchOperator.cpp)
