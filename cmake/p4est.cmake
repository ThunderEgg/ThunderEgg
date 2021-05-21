# provides imported target P4EST::P4EST
include(ExternalProject)

set(p4est_external true CACHE BOOL "build p4est library" FORCE)

# --- libsc externalProject
# this keeps libsc scope totally separate from p4est, which avoids 
# tricky to diagnose behaviors

if(NOT DEFINED P4EST_ROOT)
  set(P4EST_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  set(P4EST_LIBRARIES ${P4EST_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}p4est${CMAKE_SHARED_LIBRARY_SUFFIX})
  list(APPEND P4EST_LIBRARIES ${P4EST_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(P4EST_LIBRARIES ${P4EST_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}p4est${CMAKE_STATIC_LIBRARY_SUFFIX})
  list(APPEND P4EST_LIBRARIES ${P4EST_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(P4EST_INCLUDE_DIRS ${P4EST_ROOT}/include)

ExternalProject_Add(P4EST
GIT_REPOSITORY https://github.com/cburstedde/p4est
GIT_TAG        feature-cmake
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${P4EST_ROOT} -Dmpi:BOOL=TRUE -Dopenmp:BOOL=FALSE
BUILD_BYPRODUCTS ${P4EST_LIBRARIES}
)

# --- imported target

file(MAKE_DIRECTORY ${P4EST_INCLUDE_DIRS})
# avoid race condition

find_package(ZLIB REQUIRED)

# this GLOBAL is required to be visible via other 
# project's FetchContent of this project
add_library(P4EST::P4EST INTERFACE IMPORTED GLOBAL)
target_include_directories(P4EST::P4EST INTERFACE "${P4EST_INCLUDE_DIRS}")
target_link_libraries(P4EST::P4EST INTERFACE "${P4EST_LIBRARIES}" "${ZLIB_LIBRARIES}")

add_dependencies(P4EST::P4EST P4EST)
