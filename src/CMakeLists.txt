# -- add library

add_library(ThunderEgg)

set_target_properties(
  ThunderEgg PROPERTIES
    SOVERSION ${CMAKE_PROJECT_VERSION_MAJOR}
    VERSION   ${CMAKE_PROJECT_VERSION}
)
add_subdirectory(ThunderEgg)

# -- generate the single header ThunderEgg.h

set(ThunderEgg_HDRS_String "")
foreach(HDR ${ALL_ThunderEgg_HDRS})
  set(ThunderEgg_HDRS_String "${ThunderEgg_HDRS_String}\n#include<${HDR}>")
endforeach(HDR)
configure_file("ThunderEgg.h.in" "ThunderEgg.h")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ThunderEgg.h"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT   ThunderEgg_Development
)


# -- set include directories

target_include_directories(
  ThunderEgg PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# -- link libraries

if(TARGET P4EST::P4EST)
  target_link_libraries(ThunderEgg PUBLIC P4EST::P4EST SC::SC)
endif(TARGET P4EST::P4EST)

if(TARGET FFTW::FFTW)
  target_link_libraries(ThunderEgg PUBLIC FFTW::FFTW)
endif(TARGET FFTW::FFTW)

if(TARGET PETSC::PETSC)
  target_link_libraries(ThunderEgg PUBLIC PETSC::PETSC)
endif(TARGET PETSC::PETSC)

if(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)
  target_link_libraries(ThunderEgg PUBLIC LAPACK::LAPACK)
  target_link_libraries(ThunderEgg PUBLIC BLAS::BLAS)
endif(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)

target_link_libraries(ThunderEgg PUBLIC MPI::MPI_CXX)

# -- imported target, for use from FetchContent

add_library(ThunderEgg::ThunderEgg INTERFACE IMPORTED GLOBAL)
target_link_libraries(ThunderEgg::ThunderEgg INTERFACE ThunderEgg)