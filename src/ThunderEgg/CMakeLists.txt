list(APPEND ThunderEgg_HDRS BiLinearGhostFiller.h)
target_sources(ThunderEgg PRIVATE BiLinearGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS BiQuadraticGhostFiller.h)
target_sources(ThunderEgg PRIVATE BiQuadraticGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS BreakdownError.h)

list(APPEND ThunderEgg_HDRS BufferReader.h)

list(APPEND ThunderEgg_HDRS BufferWriter.h)

list(APPEND ThunderEgg_HDRS CoarseNbrInfo.h)

list(APPEND ThunderEgg_HDRS DivergenceError.h)

list(APPEND ThunderEgg_HDRS DimensionalArray.h)

list(APPEND ThunderEgg_HDRS Domain.h)
target_sources(ThunderEgg PRIVATE Domain.cpp)

list(APPEND ThunderEgg_HDRS DomainGenerator.h)

list(APPEND ThunderEgg_HDRS DomainTools.h)

list(APPEND ThunderEgg_HDRS Face.h)
target_sources(ThunderEgg PRIVATE Face.cpp)

list(APPEND ThunderEgg_HDRS FineNbrInfo.h)

list(APPEND ThunderEgg_HDRS GhostFiller.h)

list(APPEND ThunderEgg_HDRS Loops.h)

list(APPEND ThunderEgg_HDRS LocalData.h)
target_sources(ThunderEgg PRIVATE LocalData.cpp)

list(APPEND ThunderEgg_HDRS LocalDataManager.h)

list(APPEND ThunderEgg_HDRS MPIGhostFiller.h)
target_sources(ThunderEgg PRIVATE MPIGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS NormalNbrInfo.h)

list(APPEND ThunderEgg_HDRS NbrType.h)

list(APPEND ThunderEgg_HDRS NbrInfo.h)

list(APPEND ThunderEgg_HDRS NbrInfoBase.h)

list(APPEND ThunderEgg_HDRS Operator.h)

list(APPEND ThunderEgg_HDRS Orthant.h)
target_sources(ThunderEgg PRIVATE Orthant.cpp)

if(TARGET P4EST::P4EST)

  list(APPEND ThunderEgg_HDRS P4estDomainGenerator.h)
  target_sources(ThunderEgg PRIVATE P4estDomainGenerator.cpp)

  list(APPEND ThunderEgg_HDRS P8estDomainGenerator.h)
  target_sources(ThunderEgg PRIVATE P8estDomainGenerator.cpp)

endif(TARGET P4EST::P4EST)

list(APPEND ThunderEgg_HDRS PatchInfo.h)
target_sources(ThunderEgg PRIVATE PatchInfo.cpp)

list(APPEND ThunderEgg_HDRS PatchOperator.h)
target_sources(ThunderEgg PRIVATE PatchOperator.cpp)

list(APPEND ThunderEgg_HDRS PatchSolver.h)

list(APPEND ThunderEgg_HDRS RuntimeError.h)

list(APPEND ThunderEgg_HDRS Serializable.h)

list(APPEND ThunderEgg_HDRS Timer.h)
target_sources(ThunderEgg PRIVATE Timer.cpp)

list(APPEND ThunderEgg_HDRS TriLinearGhostFiller.h)
target_sources(ThunderEgg PRIVATE TriLinearGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS TypeDefs.h)

list(APPEND ThunderEgg_HDRS ValVector.h)
target_sources(ThunderEgg PRIVATE ValVector.cpp)

list(APPEND ThunderEgg_HDRS ValVectorGenerator.h)
target_sources(ThunderEgg PRIVATE ValVectorGenerator.cpp)

list(APPEND ThunderEgg_HDRS Vector.h)
target_sources(ThunderEgg PRIVATE Vector.cpp)

list(APPEND ThunderEgg_HDRS VectorGenerator.h)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg
)

# -- subdirectories
add_subdirectory(Experimental)
add_subdirectory(GMG)
add_subdirectory(PETSc)
add_subdirectory(Poisson)
add_subdirectory(Schur)
add_subdirectory(tpl)
add_subdirectory(VarPoisson)
add_subdirectory(Iterative)

# -- add list of headers to all

list(TRANSFORM ThunderEgg_HDRS PREPEND "ThunderEgg/")
list(APPEND ALL_ThunderEgg_HDRS ${ThunderEgg_HDRS})

# -- propogate to parent scope

set(ALL_ThunderEgg_HDRS ${ALL_ThunderEgg_HDRS} PARENT_SCOPE)


