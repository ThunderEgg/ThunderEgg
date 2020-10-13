if(PETSC_FOUND)
  add_executable(unit_tests_mpiany test.cpp OctTree.cpp LocalData.cpp Utils.cpp
                                   utils/DomainReader.cpp)
  target_include_directories(unit_tests_mpiany
                             PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(unit_tests_mpiany ThunderEgg)
  if(VTK_FOUND)
    target_link_libraries(unit_tests_mpiany VtkWriter)
  endif(VTK_FOUND)

  include(CTest)
  include(ParseAndAddCatchTests)

  set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)

  foreach(ARG1 RANGE 1 1)
    set(OptionalCatchTestLauncher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ARG1})
    set(AdditionalCatchParameters -a WORKING_DIRECTORY
                                  ${CMAKE_CURRENT_SOURCE_DIR})
    parseandaddcatchtests(
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
    unset(OptionalCatchTestLauncher)
    set(AdditionalCatchParameters)
  endforeach(ARG1 RANGE 3)

  function(add_mpi_test ARG1)
    file(GLOB_RECURSE TEST_MPI${ARG1}_SRC ${PROJECT_SOURCE_DIR}/test
         *_MPI${ARG1}.cpp)
    set(TEST_MPI${ARG1}_PETSC_SRC ${TEST_MPI${ARG1}_SRC})
    set(TEST_MPI${ARG1}_FFTW_SRC ${TEST_MPI${ARG1}_SRC})

    list(FILTER TEST_MPI${ARG1}_SRC EXCLUDE REGEX ".*_PETSC_.*")
    if(NOT FFTW_FOUND)
      list(FILTER TEST_MPI${ARG1}_SRC EXCLUDE REGEX ".*_FFTW_.*")
    endif(NOT FFTW_FOUND)

    list(FILTER TEST_MPI${ARG1}_PETSC_SRC INCLUDE REGEX ".*_PETSC_.*")

    set(OptionalCatchTestLauncher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ARG1})
    set(AdditionalCatchParameters -a WORKING_DIRECTORY
                                  ${CMAKE_CURRENT_SOURCE_DIR})

    if(TEST_MPI${ARG1}_SRC)
      add_executable(unit_tests_mpi${ARG1} ${TEST_MPI${ARG1}_SRC} test.cpp
                                           utils/DomainReader.cpp)
      target_link_libraries(unit_tests_mpi${ARG1} ThunderEgg)
      target_include_directories(unit_tests_mpi${ARG1}
                                 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
      parseandaddcatchtests(unit_tests_mpi${ARG1})
    endif()

    if(TEST_MPI${ARG1}_PETSC_SRC AND PETSC_FOUND)
      add_executable(
        unit_tests_petsc_mpi${ARG1} ${TEST_MPI${ARG1}_PETSC_SRC} test_PETSC.cpp
                                    utils/DomainReader.cpp)
      target_link_libraries(unit_tests_petsc_mpi${ARG1} ThunderEgg)
      target_include_directories(unit_tests_petsc_mpi${ARG1}
                                 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
      parseandaddcatchtests(unit_tests_petsc_mpi${ARG1})
    endif()

    unset(OptionalCatchTestLauncher)
    set(AdditionalCatchParameters)

  endfunction(add_mpi_test)

  add_mpi_test(1)
  add_mpi_test(2)
  add_mpi_test(3)

endif(PETSC_FOUND)
