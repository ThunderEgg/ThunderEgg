list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/BiCGStabSolver.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/BilinearInterpolator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/BilinearInterpolator.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/IfaceInterp.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/IfaceSet.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/IfaceType.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurHelper.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/SchurHelper.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurInfo.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/SchurInfo.cpp)

if(PETSC_FOUND)
  list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/PolyChebPrec.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/PolyChebPrec.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurMatrixHelper.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/SchurMatrixHelper.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurMatrixHelper2d.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/SchurMatrixHelper2d.cpp)
endif(PETSC_FOUND)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurWrapOp.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/TriLinInterp.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/TriLinInterp.cpp)
