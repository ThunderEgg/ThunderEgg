# -- clear header list form parent scope

unset(ThunderEgg_HDRS)

# -- add sources

list(APPEND ThunderEgg_HDRS Cycle.h)
target_sources(ThunderEgg PRIVATE Cycle.cpp)

list(APPEND ThunderEgg_HDRS CycleBuilder.h)
target_sources(ThunderEgg PRIVATE CycleBuilder.cpp)

list(APPEND ThunderEgg_HDRS CycleOpts.h)

list(APPEND ThunderEgg_HDRS DirectInterpolator.h)
target_sources(ThunderEgg PRIVATE DirectInterpolator.cpp)

list(APPEND ThunderEgg_HDRS InterLevelComm.h)
target_sources(ThunderEgg PRIVATE InterLevelComm.cpp)

list(APPEND ThunderEgg_HDRS Interpolator.h)

list(APPEND ThunderEgg_HDRS Level.h)

list(APPEND ThunderEgg_HDRS LinearRestrictor.h)
target_sources(ThunderEgg PRIVATE LinearRestrictor.cpp)

list(APPEND ThunderEgg_HDRS MPIInterpolator.h)
target_sources(ThunderEgg PRIVATE MPIInterpolator.cpp)

list(APPEND ThunderEgg_HDRS MPIRestrictor.h)
target_sources(ThunderEgg PRIVATE MPIRestrictor.cpp)

list(APPEND ThunderEgg_HDRS Restrictor.h)

list(APPEND ThunderEgg_HDRS Smoother.h)

list(APPEND ThunderEgg_HDRS VCycle.h)
target_sources(ThunderEgg PRIVATE VCycle.cpp)

list(APPEND ThunderEgg_HDRS WCycle.h)
target_sources(ThunderEgg PRIVATE WCycle.cpp)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg/GMG
)

# -- add list of headers to all

list(TRANSFORM ThunderEgg_HDRS PREPEND "ThunderEgg/GMG/")
list(APPEND ALL_ThunderEgg_HDRS ${ThunderEgg_HDRS})

# -- propogate to parent scope

set(ALL_ThunderEgg_HDRS ${ALL_ThunderEgg_HDRS} PARENT_SCOPE)