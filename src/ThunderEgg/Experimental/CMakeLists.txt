# -- clear header list form parent scope

unset(ThunderEgg_HDRS)

# -- add sources

list(APPEND ThunderEgg_HDRS OctNode.h)

list(APPEND ThunderEgg_HDRS OctTree.h)

if(TARGET PETSC::PETSC)
    list(APPEND ThunderEgg_HDRS PBMatrix.h)
    target_sources(ThunderEgg PRIVATE PBMatrix.cpp)
endif(TARGET PETSC::PETSC)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg/Experimental
)

# -- add list of headers to all

list(TRANSFORM ThunderEgg_HDRS PREPEND "ThunderEgg/Experimental/")
list(APPEND ALL_ThunderEgg_HDRS ${ThunderEgg_HDRS})

# -- propogate to parent scope

set(ALL_ThunderEgg_HDRS ${ALL_ThunderEgg_HDRS} PARENT_SCOPE)