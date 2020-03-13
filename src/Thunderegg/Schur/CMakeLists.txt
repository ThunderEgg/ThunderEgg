list(APPEND Thunderegg_HDRS Thunderegg/Schur/BiCGStabSolver.h)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/BilinearInterpolator.h)
list(APPEND Thunderegg_SRCS Thunderegg/Schur/BilinearInterpolator.cpp)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/IfaceInterp.h)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/IfaceSet.h)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/IfaceType.h)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/PatchOperator.h)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/PatchSolver.h)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/SchurHelper.h)
list(APPEND Thunderegg_SRCS Thunderegg/Schur/SchurHelper.cpp)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/SchurInfo.h)
list(APPEND Thunderegg_SRCS Thunderegg/Schur/SchurInfo.cpp)

if(PETSC_FOUND)
    list(APPEND Thunderegg_HDRS Thunderegg/Schur/PolyChebPrec.h)
    list(APPEND Thunderegg_SRCS Thunderegg/Schur/PolyChebPrec.cpp)

    list(APPEND Thunderegg_HDRS Thunderegg/Schur/SchurMatrixHelper.h)
    list(APPEND Thunderegg_SRCS Thunderegg/Schur/SchurMatrixHelper.cpp)

    list(APPEND Thunderegg_HDRS Thunderegg/Schur/SchurMatrixHelper2d.h)
    list(APPEND Thunderegg_SRCS Thunderegg/Schur/SchurMatrixHelper2d.cpp)
endif(PETSC_FOUND)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/SchurWrapOp.h)

list(APPEND Thunderegg_HDRS Thunderegg/Schur/TriLinInterp.h)
list(APPEND Thunderegg_SRCS Thunderegg/Schur/TriLinInterp.cpp)