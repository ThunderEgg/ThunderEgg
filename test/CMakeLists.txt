#cmake_policy(SET CMP0110 NEW)
include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2
  GIT_TAG        v3.0.0-preview3
)

FetchContent_MakeAvailable(catch2)


include(CTest)
include(ParseAndAddCatchTests)

add_executable(unit_tests_mpi1)
add_executable(unit_tests_mpi2)
add_executable(unit_tests_mpi3)

add_library(test_shared test.cpp)

# -- add subdirectories

add_subdirectory(GMG)
add_subdirectory(Iterative)
add_subdirectory(PETSc)
add_subdirectory(Poisson)
add_subdirectory(Schur)
add_subdirectory(utils)
add_subdirectory(VarPoisson)

# -- add sources

target_sources(unit_tests_mpi1 PRIVATE BiLinearGhostFiller_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE BiQuadraticGhostFiller_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE CoarseNbrInfo_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE ComponentView_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE ConstView_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Corner_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Domain_MPI1.cpp)
target_sources(unit_tests_mpi2 PRIVATE Domain_MPI2.cpp)

target_sources(unit_tests_mpi1 PRIVATE DomainTools_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Edge_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Face_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE FineNbrInfo_MPI1.cpp)


target_sources(unit_tests_mpi1 PRIVATE MPIGhostFiller_MPI1.cpp)
target_sources(unit_tests_mpi2 PRIVATE MPIGhostFiller_MPI2.cpp)
target_sources(unit_tests_mpi3 PRIVATE MPIGhostFiller_MPI3.cpp)

target_sources(unit_tests_mpi1 PRIVATE NbrType_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE NormalNbrInfo_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Orthant_MPI1.cpp)

if(TARGET P4EST::P4EST)

  target_sources(unit_tests_mpi1 PRIVATE P4estDomainGenerator_MPI1.cpp)
  target_sources(unit_tests_mpi2 PRIVATE P4estDomainGenerator_MPI2.cpp)

  target_sources(test_shared PRIVATE P8estDomainGenerator_Check2x2x2DomainNbrs_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check2x2x2RefinedBSWDomainNbrs_SHARED.cpp)

  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainCornerHasNbrs_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainCornerNbrIds_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainCornerNbrRanks_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainEdgeHasNbrs_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainEdgeNbrIds_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainEdgeNbrRanks_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainSideHasNbrs_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainSideNbrIds_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4DomainSideNbrRanks_SHARED.cpp)

  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainCornerHasNbrs_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainCornerNbrIds_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainCornerNbrOrths_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainCornerNbrRanks_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainEdgeHasNbrs_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainEdgeNbrIds_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainEdgeNbrOrths_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainEdgeNbrRanks_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainSideHasNbrs_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainSideNbrIds_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainSideNbrOrths_SHARED.cpp)
  target_sources(test_shared PRIVATE P8estDomainGenerator_Check4x4x4RefinedBSWDomainSideNbrRanks_SHARED.cpp)

  target_sources(unit_tests_mpi1 PRIVATE P8estDomainGenerator_MPI1.cpp)
  target_sources(unit_tests_mpi2 PRIVATE P8estDomainGenerator_MPI2.cpp)

  target_sources(test_shared PRIVATE P8estDomainGenerator_SHARED.cpp)
  
endif(TARGET P4EST::P4EST)

target_sources(unit_tests_mpi1 PRIVATE PatchInfo_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE PatchArray_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE PatchOperator_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE PatchSolver_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Side_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Timer_MPI1.cpp)
target_sources(unit_tests_mpi2 PRIVATE Timer_MPI2.cpp)

target_sources(unit_tests_mpi1 PRIVATE TriLinearGhostFiller_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE ValVector_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE ValVectorGenerator_MPI1.cpp)

target_sources(unit_tests_mpi1 PRIVATE Vector_MPI1.cpp)
target_sources(unit_tests_mpi2 PRIVATE Vector_MPI2.cpp)
target_sources(unit_tests_mpi3 PRIVATE Vector_MPI3.cpp)

target_sources(unit_tests_mpi1 PRIVATE View_MPI1.cpp)


# -- determine compile flags

if(TARGET P4EST::P4EST)
  list(APPEND TEST_COMPILE_FLAGS TEST_P4EST=1)
endif(TARGET P4EST::P4EST)

if(TARGET PETSC::PETSC)
  list(APPEND TEST_COMPILE_FLAGS TEST_PETSC=1)
endif(TARGET PETSC::PETSC)

list(APPEND TEST_COMPILE_FLAGS TEST_DIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}\")

target_compile_definitions(unit_tests_mpi1 PRIVATE ${TEST_COMPILE_FLAGS})
target_compile_definitions(unit_tests_mpi2 PRIVATE ${TEST_COMPILE_FLAGS})
target_compile_definitions(unit_tests_mpi3 PRIVATE ${TEST_COMPILE_FLAGS})
target_compile_definitions(test_shared PRIVATE ${TEST_COMPILE_FLAGS})

# -- link libraries

target_link_libraries(test_shared PUBLIC ThunderEgg Catch2)

target_link_libraries(unit_tests_mpi1 test_shared)
target_link_libraries(unit_tests_mpi2 test_shared)
target_link_libraries(unit_tests_mpi3 test_shared)

    #target_include_directories(unit_tests_mpi${ARG2}
    #                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)

function(add_mpi_test ARG1 ARG2)
  set(OptionalCatchTestLauncher ${MPIEXEC} "--oversubscribe" ${MPIEXEC_NUMPROC_FLAG}  ${ARG2})

  parseandaddcatchtests(${ARG1})

  unset(OptionalCatchTestLauncher)
  set(AdditionalCatchParameters)

endfunction(add_mpi_test)

add_mpi_test(unit_tests_mpi1 1)
add_mpi_test(unit_tests_mpi2 2)
add_mpi_test(unit_tests_mpi3 3)