# -- clear header list form parent scope

unset(ThunderEgg_HDRS)

# -- add sources

if(TARGET PETSc::PETSc)

  list(APPEND ThunderEgg_HDRS KSPSolver.h)
  list(APPEND ThunderEgg_HDRS MatWrapper.h)
  list(APPEND ThunderEgg_HDRS MatShellCreator.h)
  list(APPEND ThunderEgg_HDRS PCShellCreator.h)

endif(TARGET PETSc::PETSc)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg/PETSc
  COMPONENT   ThunderEgg_Development
)

# -- add list of headers to all

list(TRANSFORM ThunderEgg_HDRS PREPEND "ThunderEgg/PETSc/")
list(APPEND ALL_ThunderEgg_HDRS ${ThunderEgg_HDRS})

# -- propogate to parent scope

set(ALL_ThunderEgg_HDRS ${ALL_ThunderEgg_HDRS} PARENT_SCOPE)