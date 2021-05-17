# provides imported target SC::SC
include(ExternalProject)

set(sc_external true CACHE BOOL "build sc library" FORCE)

# --- libsc externalProject
# this keeps libsc scope totally separate from p4est, which avoids 
# tricky to diagnose behaviors

if(NOT DEFINED SC_ROOT)
  set(SC_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  set(SC_LIBRARIES ${SC_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(SC_LIBRARIES ${SC_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(SC_INCLUDE_DIRS ${SC_ROOT}/include)

find_package(ZLIB REQUIRED)

# --- imported target

file(MAKE_DIRECTORY ${SC_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other 
# project's FetchContent of this project
add_library(SC::SC INTERFACE IMPORTED GLOBAL)
target_include_directories(SC::SC INTERFACE "${SC_INCLUDE_DIRS}")
target_link_libraries(SC::SC INTERFACE "${SC_LIBRARIES}")
target_link_libraries(SC::SC INTERFACE "${ZLIB_LIBRARIES}")

add_dependencies(SC::SC SC)
