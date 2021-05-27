# -- determine sources and public headers first

list(APPEND ThunderEgg_HDRS ThunderEgg/BiLinearGhostFiller.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/BiLinearGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/BiQuadraticGhostFiller.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/BiQuadraticGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/BreakdownError.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/BufferReader.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/BufferWriter.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/CoarseNbrInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/DivergenceError.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/DimensionalArray.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Domain.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Domain.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/DomainGenerator.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/DomainTools.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Face.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Face.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/FineNbrInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/GhostFiller.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Loops.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/LocalData.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/LocalData.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/MPIGhostFiller.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/MPIGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/NormalNbrInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/NbrType.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/NbrInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/NbrInfoBase.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Operator.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Orthant.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Orthant.cpp)

if(TARGET P4EST::P4EST)

  list(APPEND ThunderEgg_HDRS ThunderEgg/P4estDomainGenerator.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/P4estDomainGenerator.cpp)

  list(APPEND ThunderEgg_HDRS ThunderEgg/P8estDomainGenerator.h)
  list(APPEND ThunderEgg_SRCS ThunderEgg/P8estDomainGenerator.cpp)

endif(TARGET P4EST::P4EST)

list(APPEND ThunderEgg_HDRS ThunderEgg/PatchInfo.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/PatchInfo.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/PatchOperator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/PatchOperator.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/PatchSolver.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/RuntimeError.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Serializable.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Timer.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Timer.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/TriLinearGhostFiller.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/TriLinearGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/TypeDefs.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/ValVector.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/ValVector.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/ValVectorGenerator.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/ValVectorGenerator.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/Vector.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Vector.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/VectorGenerator.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/tpl/json.hpp)

include(ThunderEgg/Experimental/CMakeLists.txt)
include(ThunderEgg/GMG/CMakeLists.txt)
include(ThunderEgg/PETSc/CMakeLists.txt)
include(ThunderEgg/Poisson/CMakeLists.txt)
include(ThunderEgg/Schur/CMakeLists.txt)
include(ThunderEgg/VarPoisson/CMakeLists.txt)
include(ThunderEgg/Iterative/CMakeLists.txt)

add_library(ThunderEgg ${ThunderEgg_SRCS})

# -- generate the single header ThunderEgg.h

set(ThunderEgg_HDRS_String "")
foreach(HDR ${ThunderEgg_HDRS})
  set(ThunderEgg_HDRS_String "${ThunderEgg_HDRS_String}\n#include<${HDR}>")
endforeach(HDR)
configure_file("ThunderEgg.h.in" "ThunderEgg.h")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ThunderEgg.h"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


# -- set include directories

target_include_directories(
  ThunderEgg INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# -- link libraries

if(TARGET P4EST::P4EST)
  target_link_libraries(ThunderEgg PUBLIC P4EST::P4EST)
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

# -- install

install(
  TARGETS ThunderEgg
  EXPORT ThunderEggConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
}

install(EXPORT ThunderEggConfig DESTINATION share/ThunderEgg/cmake)

# -- generate a pkgconfig file

include(PkgConfigHandler)
configure_pc_file("ThunderEgg" "ThunderEgg.pc.in" "ThunderEgg.pc")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ThunderEgg.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# -- imported target, for use from FetchContent

add_library(ThunderEgg::ThunderEgg INTERFACE IMPORTED GLOBAL)
target_link_libraries(ThunderEgg::ThunderEgg INTERFACE ThunderEgg)