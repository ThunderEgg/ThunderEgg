cmake_minimum_required(VERSION 3.18)
project(ThunderEgg 
        VERSION 1.0.2
        LANGUAGES C CXX)
enable_testing()

# set(CMAKE_CXX_EXTENSIONS OFF)
include(GNUInstallDirs)

include(cmake/options.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

find_package(MPI REQUIRED COMPONENTS C CXX)

set(THUNDEREGG_ENABLED_COMPONENTS "")

if(petsc)

  find_package(PETSc)

  if(petsc_required AND NOT PETSc_FOUND)
    message(FATAL_ERROR "PETSc was not found")
  endif()

  if(PETSc_FOUND)
    set(THUNDEREGG_PETSC_ENABLED TRUE)
    list(APPEND THUNDEREGG_ENABLED_COMPONENTS "PETSC")
  endif()

endif()

if(fftw)

  find_package(FFTW)

  if(fftw_required AND NOT FFTW_FOUND)
    message(FATAL_ERROR "FFTW was not found")
  endif()

  if(FFTW_FOUND)
    set(THUNDEREGG_FFTW_ENABLED TRUE)
    list(APPEND THUNDEREGG_ENABLED_COMPONENTS "FFTW")
  endif()
endif()

if(lapack)

  find_package(LAPACK)
  find_package(BLAS)

  if(lapack_required AND NOT LAPACK_FOUND)
    message(FATAL_ERROR "LAPACK was not found")
  endif()
  if(lapack_required AND NOT BLAS_FOUND)
    message(FATAL_ERROR "BLAS was not found")
  endif()

  if(BLAS_FOUND AND LAPACK_FOUND)
    set(THUNDEREGG_LAPACK_ENABLED TRUE)
    list(APPEND THUNDEREGG_ENABLED_COMPONENTS "LAPACK")
  endif()
endif()


if(p4est)

  find_package(P4EST)
  find_package(SC)

  if(NOT (P4EST_FOUND AND SC_FOUND) AND p4est_external)
    include(cmake/p4est.cmake)
  endif()

  if(p4est_required AND NOT TARGET P4EST::P4EST)
    message(FATAL_ERROR "p4est was not found")
  endif()
  if(P4EST_FOUND AND NOT SC_FOUND)
    message(FATAL_ERROR "SC was not found")
  endif()

  if(TARGET P4EST::P4EST)
    set(THUNDEREGG_P4EST_ENABLED TRUE)
    list(APPEND THUNDEREGG_ENABLED_COMPONENTS "P4EST")
  endif()
endif()

add_subdirectory(src)


include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)

include(cmake/documentation.cmake)

include(CMakePackageConfigHelpers)


configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION lib
  NO_CHECK_REQUIRED_COMPONENTS_MACRO  
)

  install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindFFTW.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindPETSc.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindP4EST.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindSC.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/CorrectWindowsPaths.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/ResolveCompilerPaths.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindPackageMultipass.cmake
  DESTINATION lib/cmake/Modules)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  COMPATIBILITY SameMinorVersion)

install(TARGETS ThunderEgg 
  EXPORT                        ThunderEggTargets
  RUNTIME DESTINATION           ${CMAKE_INSTALL_BINDIR}
          COMPONENT             ThunderEgg_RunTime
  LIBRARY DESTINATION           ${CMAKE_INSTALL_LIBDIR}
          COMPONENT             ThunderEgg_RunTime
          NAMELINK_COMPONENT    ThunderEgg_Development
  ARCHIVE DESTINATION           ${CMAKE_INSTALL_LIBDIR}
          COMPONENT             ThunderEgg_Development
)

install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  COMPONENT   ThunderEgg_Development
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  COMPONENT   ThunderEgg_Development
)

export(EXPORT ${PROJECT_NAME}Targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
)

# -- generate a pkgconfig file
include(cmake/pkgconfig.cmake)
