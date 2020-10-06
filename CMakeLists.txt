cmake_minimum_required(VERSION 2.8.8)
project(ThunderEgg)
enable_language(CXX)
enable_language(C)
enable_testing()

# set(CMAKE_CXX_EXTENSIONS OFF)
include(GNUInstallDirs)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3") set
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wall -g -O3")
set(CMAKE_CXX_STANDARD 11)

# find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message("-- Using ccache.")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

find_package(MPI REQUIRED)
find_package(PETSc)
find_package(FFTW)
find_package(Zoltan)
find_package(p4est)
find_package(BLAS)
find_package(LAPACK)
if(p4est_FOUND)
  find_package(sc REQUIRED)
endif()

add_subdirectory(tpl)
add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(test)

find_package(Doxygen)
option(BUILD_DOCUMENTATION
       "Create and install the HTML based API documentation (requires Doxygen)"
       ${DOXYGEN_FOUND})

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
