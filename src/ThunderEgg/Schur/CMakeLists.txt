list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/CoarseIfaceInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/IfaceInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/IfaceType.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/Interface.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/Interface.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/NormalIfaceInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/InterfaceDomain.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/InterfaceDomain.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/PatchIfaceInfo.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/PatchIfaceInfo.cpp)

if(PETSC_FOUND)
  list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/PolyChebPrec.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/PolyChebPrec.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurMatrixHelper.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/SchurMatrixHelper.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurMatrixHelper2d.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/SchurMatrixHelper2d.cpp)
endif(PETSC_FOUND)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/SchurWrapOp.h)
