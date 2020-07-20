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
  include(Catch)

  foreach(ARG1 RANGE 1 1)
    catch_discover_tests(
      unit_tests_mpiany
      WORKING_DIRECTORY
      ${PROJECT_SOURCE_DIR}/test
      THUNDEREGG_TEST_EXECUTOR
      ${MPIEXEC}
      ${MPIEXEC_NUMPROC_FLAG}
      ${ARG1}
      THUNDEREGG_SINGLE_EXECUTOR
      ${MPIEXEC}
      ${MPIEXEC_NUMPROC_FLAG}
      1
      -q)
  endforeach(ARG1 RANGE 3)

  function(add_mpi_test ARG1)
    file(GLOB_RECURSE TEST_MPI${ARG1}_SRC ${PROJECT_SOURCE_DIR}/test
         *_MPI${ARG1}.cpp)
    set(TEST_MPI${ARG1}_PETSC_SRC ${TEST_MPI${ARG1}_SRC})

    list(FILTER TEST_MPI${ARG1}_SRC EXCLUDE REGEX ".*_PETSC_.*")
    list(FILTER TEST_MPI${ARG1}_PETSC_SRC INCLUDE REGEX ".*_PETSC_.*")

    if(TEST_MPI${ARG1}_SRC)
      add_executable(unit_tests_mpi${ARG1} ${TEST_MPI${ARG1}_SRC} test.cpp
                                           utils/DomainReader.cpp)
      target_link_libraries(unit_tests_mpi${ARG1} Thunderegg)
      target_include_directories(unit_tests_mpi${ARG1}
                                 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
      catch_discover_tests(
        unit_tests_mpi${ARG1}
        WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}/test
        THUNDEREGG_TEST_EXECUTOR
        ${MPIEXEC}
        ${MPIEXEC_NUMPROC_FLAG}
        ${ARG1}
        THUNDEREGG_SINGLE_EXECUTOR
        ${MPIEXEC}
        ${MPIEXEC_NUMPROC_FLAG}
        1
        -q)
    endif()
    if(TEST_MPI${ARG1}_PETSC_SRC AND PETSC_FOUND)
      add_executable(
        unit_tests_petsc_mpi${ARG1} ${TEST_MPI${ARG1}_PETSC_SRC} test_PETSC.cpp
                                    utils/DomainReader.cpp)
      target_link_libraries(unit_tests_petsc_mpi${ARG1} Thunderegg)
      target_include_directories(unit_tests_petsc_mpi${ARG1}
                                 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
      catch_discover_tests(
        unit_tests_petsc_mpi${ARG1}
        WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}/test
        EXTRA_ARGS
        -x
        10
        THUNDEREGG_TEST_EXECUTOR
        ${MPIEXEC}
        ${MPIEXEC_NUMPROC_FLAG}
        ${ARG1}
        THUNDEREGG_SINGLE_EXECUTOR
        ${MPIEXEC}
        ${MPIEXEC_NUMPROC_FLAG}
        1
        -q)
    endif()
  endfunction(add_mpi_test)

  add_mpi_test(1)
  add_mpi_test(2)
  add_mpi_test(3)

endif(PETSC_FOUND)
