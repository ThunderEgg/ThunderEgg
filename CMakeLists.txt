cmake_minimum_required(VERSION 3.18)
project(ThunderEgg 
        VERSION 1.0.1
        LANGUAGES CXX)
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

if(NOT disable_petsc)
  find_package(PETSc)
endif(NOT disable_petsc)

find_package(FFTW)
find_package(BLAS)
find_package(LAPACK)
find_package(P4EST)
if(TARGET P4EST::P4EST)
  find_package(SC REQUIRED HINTS ${P4EST_ROOT})
  enable_language(C)
  find_package(MPI COMPONENTS C REQRUIRED)
else()
  include(cmake/p4est.cmake)
endif()

add_subdirectory(src)


include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)

include(cmake/documentation.cmake)

include(CMakePackageConfigHelpers)

if(TARGET FFTW::FFTW)
  set(THUNDEREGG_FFTW_ENABLED TRUE)
endif()
if(TARGET P4EST::P4EST)
  set(THUNDEREGG_P4EST_ENABLED TRUE)
endif()
if(TARGET BLAS::BLAS AND TARGET LAPACK::LAPACK)
  set(THUNDEREGG_LAPACK_ENABLED TRUE)
endif()
if(TARGET PETSC::PETSC)
  set(THUNDEREGG_PETSC_ENABLED TRUE)
endif()

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION lib)

  install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindFFTW.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindPETSc.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/CorrectWindowsPaths.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/ResolveCompilerPaths.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindPackageMultipass.cmake
  DESTINATION lib/cmake/Modules)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  COMPATIBILITY SameMinorVersion)

install(TARGETS ThunderEgg EXPORT ThunderEggTargets)

install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION lib/cmake/${PROJECT_NAME})

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})

export(EXPORT ${PROJECT_NAME}Targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::)

# -- generate a pkgconfig file
include(cmake/pkgconfig.cmake)