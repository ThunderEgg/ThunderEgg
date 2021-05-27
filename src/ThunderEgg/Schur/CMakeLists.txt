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

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/PatchIfaceScatter.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/PatchIfaceScatter.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/PatchSolverWrapper.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Schur/PatchSolverWrapper.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/ValVectorGenerator.h)

if(TARGET PETSC::PETSC)
  list(APPEND ThunderEgg_HDRS ThunderEgg/Schur/VecWrapperGenerator.h)
endif(TARGET PETSC::PETSC)
