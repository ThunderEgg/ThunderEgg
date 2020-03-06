if(PETSC_FOUND)
    list(APPEND Thunderegg_HDRS Thunderegg/Poisson/MatrixHelper.h)
    list(APPEND Thunderegg_SRCS Thunderegg/Poisson/MatrixHelper.cpp)

    list(APPEND Thunderegg_HDRS Thunderegg/Poisson/MatrixHelper2d.h)
    list(APPEND Thunderegg_SRCS Thunderegg/Poisson/MatrixHelper2d.cpp)

    list(APPEND Thunderegg_HDRS Thunderegg/Poisson/StencilHelper.h)
    list(APPEND Thunderegg_HDRS Thunderegg/Poisson/StencilHelper2d.h)
endif(PETSC_FOUND)

list(APPEND Thunderegg_HDRS Thunderegg/Poisson/Schur/DftPatchSolver.h)
list(APPEND Thunderegg_SRCS Thunderegg/Poisson/Schur/DftPatchSolver.cpp)

if(FFTW_FOUND)

    list(APPEND Thunderegg_HDRS Thunderegg/Poisson/Schur/FftwPatchSolver.h)
    list(APPEND Thunderegg_SRCS Thunderegg/Poisson/Schur/FftwPatchSolver.cpp)

endif(FFTW_FOUND)

list(APPEND Thunderegg_HDRS Thunderegg/Poisson/Schur/SevenPtPatchOperator.h)
list(APPEND Thunderegg_SRCS Thunderegg/Poisson/Schur/SevenPtPatchOperator.cpp)

list(APPEND Thunderegg_HDRS Thunderegg/Poisson/Schur/StarPatchOperator.h)
list(APPEND Thunderegg_SRCS Thunderegg/Poisson/Schur/StarPatchOperator.cpp)