cmake_minimum_required(VERSION 2.8.8)
project(Thunderegg)
enable_language(CXX)
enable_language(C)
enable_testing()

# set(CMAKE_CXX_EXTENSIONS OFF)
include(GNUInstallDirs)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3") set
# (CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -DDD_DEBUG")
# set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wall
# -g -O3 -DDD_DEBUG")
set(CMAKE_CXX_STANDARD 11)

# find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message("-- Using ccache.")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

find_package(MPI REQUIRED)

add_subdirectory(tpl)
add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(test)
