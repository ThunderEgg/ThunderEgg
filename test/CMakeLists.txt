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

set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)

file(GLOB_RECURSE TEST_SHARED_SRC ${PROJECT_SOURCE_DIR}/test *_SHARED.cpp)
add_library(test_shared ${TEST_SHARED_SRC})
 target_link_libraries(test_shared ThunderEgg Catch2)

function(add_mpi_test ARG1)
  file(GLOB_RECURSE TEST_MPI${ARG1}_SRC ${PROJECT_SOURCE_DIR}/test
       *_MPI${ARG1}.cpp)
  set(TEST_MPI${ARG1}_PETSC_SRC ${TEST_MPI${ARG1}_SRC})
  set(TEST_MPI${ARG1}_FFTW_SRC ${TEST_MPI${ARG1}_SRC})

  list(FILTER TEST_MPI${ARG1}_SRC EXCLUDE REGEX ".*_PETSC_.*")
  if(NOT TARGET FFTW::FFTW)
    list(FILTER TEST_MPI${ARG1}_SRC EXCLUDE REGEX ".*_FFTW_.*")
  endif(NOT TARGET FFTW::FFTW)
  if(NOT TARGET P4EST::P4EST)
    list(FILTER TEST_MPI${ARG1}_SRC EXCLUDE REGEX ".*_P4EST_.*")
  endif(NOT TARGET P4EST::P4EST)

  list(FILTER TEST_MPI${ARG1}_PETSC_SRC INCLUDE REGEX ".*_PETSC_.*")

  set(OptionalCatchTestLauncher ${MPIEXEC} "--oversubscribe" ${MPIEXEC_NUMPROC_FLAG}  ${ARG1})
  set(AdditionalCatchParameters -a WORKING_DIRECTORY
                                ${CMAKE_CURRENT_SOURCE_DIR})

  if(TEST_MPI${ARG1}_SRC)
    add_executable(unit_tests_mpi${ARG1} ${TEST_MPI${ARG1}_SRC} test.cpp)
    target_link_libraries(unit_tests_mpi${ARG1} test_shared)
    target_include_directories(unit_tests_mpi${ARG1}
                               PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    if(TARGET P4EST::P4EST)
      set_source_files_properties(test.cpp PROPERTIES COMPILE_FLAGS "-DTEST_P4EST")
    ENDIF(TARGET P4EST::P4EST)
    parseandaddcatchtests(unit_tests_mpi${ARG1})
  endif()

  if(TEST_MPI${ARG1}_PETSC_SRC AND TARGET PETSC::PETSC)
    add_executable(
      unit_tests_petsc_mpi${ARG1} ${TEST_MPI${ARG1}_PETSC_SRC} test_PETSC.cpp)
    
    target_link_libraries(unit_tests_petsc_mpi${ARG1} test_shared)
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

