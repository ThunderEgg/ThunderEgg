cmake_minimum_required(VERSION 3.14)
project(ThunderEgg LANGUAGES C CXX)
enable_testing()

# set(CMAKE_CXX_EXTENSIONS OFF)
include(GNUInstallDirs)

include(cmake/options.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message("-- Using ccache.")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)
find_package(MPI REQUIRED)
find_package(PETSc)
find_package(FFTW)
find_package(BLAS)
find_package(LAPACK)

if(NOT TARGET P4EST::P4EST)
  if(NOT p4est_external)
    find_package(P4EST HINTS ${P4EST_DIR})
  endif(NOT p4est_external)
  if(p4est_external OR NOT TARGET P4EST::P4EST)
    include(FetchContent)
    FetchContent_Declare(
      p4est
      GIT_REPOSITORY https://github.com/cburstedde/p4est
      GIT_TAG        feature-cmake
    )
    FetchContent_MakeAvailable(p4est)
  endif()
endif(NOT TARGET P4EST::P4EST)
if(TARGET P4EST::P4EST AND NOT TARGET SC::SC)
  find_package(SC PATHS ${P4EST_DIR} ${SC_DIR} REQUIRED)
  find_package(ZLIB REQUIRED)
endif(TARGET P4EST::P4EST AND NOT TARGET SC::SC)

add_subdirectory(tpl)
add_subdirectory(src)


if(examples)
  add_subdirectory(apps)
endif(examples)
if(unit_tests)
  add_subdirectory(test)
endif(unit_tests)

find_package(Doxygen)
if(Doxygen_FOUND)
  option(BUILD_DOCUMENTATION
       "Create and install the HTML based API documentation (requires Doxygen)"
       ${DOXYGEN_FOUND})
endif(Doxygen_FOUND)

if(BUILD_DOCUMENTATION)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()

  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

  add_custom_target(
    doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif()
