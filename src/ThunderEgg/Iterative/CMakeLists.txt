# -- clear header list form parent scope

unset(ThunderEgg_HDRS)

# -- add sources

list(APPEND ThunderEgg_HDRS BiCGStab.h)
target_sources(ThunderEgg PRIVATE BiCGStab.cpp)

list(APPEND ThunderEgg_HDRS BreakdownError.h)

list(APPEND ThunderEgg_HDRS CG.h)
target_sources(ThunderEgg PRIVATE CG.cpp)

list(APPEND ThunderEgg_HDRS PatchSolver.h)
target_sources(ThunderEgg PRIVATE PatchSolver.cpp)

list(APPEND ThunderEgg_HDRS Solver.h)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg/Iterative
)

# -- add list of headers to all

list(TRANSFORM ThunderEgg_HDRS PREPEND "ThunderEgg/Iterative/")
list(APPEND ALL_ThunderEgg_HDRS ${ThunderEgg_HDRS})

# -- propogate to parent scope

set(ALL_ThunderEgg_HDRS ${ALL_ThunderEgg_HDRS} PARENT_SCOPE)