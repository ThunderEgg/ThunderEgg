# -- clear header list form parent scope

unset(ThunderEgg_HDRS)

# -- add sources

list(APPEND ThunderEgg_HDRS CoarseIfaceInfo.h)

list(APPEND ThunderEgg_HDRS FineIfaceInfo.h)

list(APPEND ThunderEgg_HDRS IfaceInfo.h)

list(APPEND ThunderEgg_HDRS IfaceType.h)

list(APPEND ThunderEgg_HDRS Interface.h)
target_sources(ThunderEgg PRIVATE Interface.cpp)

list(APPEND ThunderEgg_HDRS InterfaceDomain.h)
target_sources(ThunderEgg PRIVATE InterfaceDomain.cpp)

list(APPEND ThunderEgg_HDRS NormalIfaceInfo.h)

list(APPEND ThunderEgg_HDRS PatchIfaceInfo.h)
target_sources(ThunderEgg PRIVATE PatchIfaceInfo.cpp)

list(APPEND ThunderEgg_HDRS PatchIfaceScatter.h)
target_sources(ThunderEgg PRIVATE PatchIfaceScatter.cpp)

list(APPEND ThunderEgg_HDRS PatchSolverWrapper.h)
target_sources(ThunderEgg PRIVATE PatchSolverWrapper.cpp)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg/Schur
)

# -- add list of headers to all

list(TRANSFORM ThunderEgg_HDRS PREPEND "ThunderEgg/Schur/")
list(APPEND ALL_ThunderEgg_HDRS ${ThunderEgg_HDRS})

# -- propogate to parent scope

set(ALL_ThunderEgg_HDRS ${ALL_ThunderEgg_HDRS} PARENT_SCOPE)