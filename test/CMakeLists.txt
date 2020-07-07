if(PETSC_FOUND)
  add_executable(
    unit_tests_mpiany
    test.cpp
    SchurInfo.cpp
    PatchInfo.cpp
    test.cpp
    OctTree.cpp
    Domain.cpp
    LocalData.cpp
    DomainTools.cpp
    SimpleGhostFiller.cpp
    Utils.cpp
    GMG/CycleFactory.cpp
    utils/DomainReader.cpp
    Timer.cpp)
  target_include_directories(unit_tests_mpiany
                             PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(unit_tests_mpiany Thunderegg)
  if(VTK_FOUND)
    target_link_libraries(unit_tests_mpiany VtkWriter)
  endif(VTK_FOUND)

  include(CTest)
  include(ParseAndAddCatchTests)
  set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)

  foreach(ARG1 RANGE 1 1)
    set(OptionalCatchTestLauncher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ARG1})
    set(AdditionalCatchParameters WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    parseandaddcatchtests(unit_tests_mpiany)
    unset(OptionalCatchTestLauncher)
    unset(AdditionalCatchParameters)

  endforeach(ARG1 RANGE 3)

  function(add_mpi_test ARG1)
    set(OptionalCatchTestLauncher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ARG1})
    set(AdditionalCatchParameters WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    file(GLOB_RECURSE TEST_MPI${ARG1}_SRC ${PROJECT_SOURCE_DIR}/test
         *_MPI${ARG1}.cpp)
    add_executable(unit_tests_mpi${ARG1} ${TEST_MPI${ARG1}_SRC} test.cpp
                                         utils/DomainReader.cpp)
    target_link_libraries(unit_tests_mpi${ARG1} Thunderegg)
    target_include_directories(unit_tests_mpi${ARG1}
                               PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    parseandaddcatchtests(unit_tests_mpi${ARG1})
    unset(OptionalCatchTestLauncher)
    unset(AdditionalCatchParameters)
  endfunction(add_mpi_test)

  add_mpi_test(1)
  add_mpi_test(2)
  add_mpi_test(3)

endif(PETSC_FOUND)
