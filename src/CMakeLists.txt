project(Thunderegg_Lib)

# determine sources first
list(APPEND Thunderegg_HDRS Thunderegg/BiCGStab.h)

list(APPEND Thunderegg_HDRS Thunderegg/BiCGStabPatchSolver.h)

list(APPEND Thunderegg_HDRS Thunderegg/BiLinearGhostFiller.h)

list(APPEND Thunderegg_HDRS Thunderegg/BufferIO.h)

list(APPEND Thunderegg_HDRS Thunderegg/Domain.h)
list(APPEND Thunderegg_SRCS Thunderegg/Domain.cpp)

list(APPEND Thunderegg_HDRS Thunderegg/DomainGenerator.h)

list(APPEND Thunderegg_HDRS Thunderegg/DomainTools.h)

list(APPEND Thunderegg_HDRS Thunderegg/DomainWrapOp.h)

list(APPEND Thunderegg_HDRS Thunderegg/GhostFiller.h)

list(APPEND Thunderegg_HDRS Thunderegg/Loops.h)

list(APPEND Thunderegg_HDRS Thunderegg/Operator.h)

if(p4est_FOUND)

  list(APPEND Thunderegg_HDRS Thunderegg/P4estDomGen.h)
  list(APPEND Thunderegg_HDRS Thunderegg/P4estDomGen.cpp)

endif(p4est_FOUND)

list(APPEND Thunderegg_HDRS Thunderegg/PatchInfo.h)
list(APPEND Thunderegg_SRCS Thunderegg/PatchInfo.cpp)

list(APPEND Thunderegg_HDRS Thunderegg/PatchOperator.h)

list(APPEND Thunderegg_HDRS Thunderegg/PatchSolver.h)

if(PETSC_FOUND)

  list(APPEND Thunderegg_HDRS Thunderegg/PetscMatOp.h)
  list(APPEND Thunderegg_HDRS Thunderegg/PetscShellCreator.h)
  list(APPEND Thunderegg_HDRS Thunderegg/PetscVector.h)
  list(APPEND Thunderegg_HDRS Thunderegg/PW.h)

endif(PETSC_FOUND)

list(APPEND Thunderegg_HDRS Thunderegg/SchurDomainOp.h)

list(APPEND Thunderegg_HDRS Thunderegg/SchwarzPrec.h)

list(APPEND Thunderegg_HDRS Thunderegg/Serializable.h)

list(APPEND Thunderegg_HDRS Thunderegg/Side.h)

list(APPEND Thunderegg_HDRS Thunderegg/SimpleGhostFiller.h)

list(APPEND Thunderegg_HDRS Thunderegg/Timer.h)

list(APPEND Thunderegg_HDRS Thunderegg/TypeDefs.h)

list(APPEND Thunderegg_HDRS Thunderegg/ValVector.h)

list(APPEND Thunderegg_HDRS Thunderegg/Vector.h)
list(APPEND Thunderegg_SRCS Thunderegg/Vector.cpp)

include(Thunderegg/Poisson/CMakeLists.txt)
include(Thunderegg/VarPoisson/CMakeLists.txt)
include(Thunderegg/GMG/CMakeLists.txt)
include(Thunderegg/Schur/CMakeLists.txt)
include(Thunderegg/Experimental/CMakeLists.txt)

add_library(Thunderegg ${Thunderegg_SRCS} ${Thunderegg_HDRS})
target_compile_features(Thunderegg PUBLIC cxx_std_11)

# foreach(HDR ${Thunderegg_HDRS}) get_filename_component(DIR ${HDR} PATH)
# file(COPY ${HDR} DESTINATION include/${DIR}) endforeach()

# set include directories
target_include_directories(
  Thunderegg PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
list(APPEND Thunderegg_Includes ${MPI_CXX_INCLUDE_PATH})
list(APPEND Thunderegg_Includes ${PETSC_INCLUDES})
list(APPEND Thunderegg_Includes ${Zoltan_INCLUDES})
if(p4est_FOUND)
  list(APPEND Thunderegg_Includes ${p4est_INCLUDES})
  list(APPEND Thunderegg_Includes ${sc_INCLUDES})
endif()
target_include_directories(Thunderegg PUBLIC ${Thunderegg_Includes})

# set link list(APPEND Thunderegg_Libs ${MPI_CXX_LIBRARIES})
list(APPEND Thunderegg_Libs ${PETSC_LIBRARIES})
if(FFTW_FOUND)
  list(APPEND Thunderegg_Libs ${FFTW_LIBRARIES})
  list(APPEND Thunderegg_Includes ${FFTW_INCLUDES})
endif(FFTW_FOUND)
list(APPEND Thunderegg_Libs ${Zoltan_LIBRARIES})
list(APPEND Thunderegg_Libs ${BLAS_LIBRARIES})
list(APPEND Thunderegg_Libs ${LAPACK_LIBRARIES})
if(p4est_FOUND)
  list(APPEND Thunderegg_Libs ${p4est_LIBRARIES})
  list(APPEND Thunderegg_Libs ${sc_LIBRARIES})
endif()
target_link_libraries(Thunderegg PUBLIC ${Thunderegg_Libs})
target_link_libraries(Thunderegg PUBLIC ${MPI_CXX_LIBRARIES})

install(
  TARGETS Thunderegg
  EXPORT ThundereggConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(
  DIRECTORY Thunderegg/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Thunderegg
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*.hpp")

install(EXPORT ThundereggConfig DESTINATION share/Thunderegg/cmake)

include(PkgConfigHandler)
configure_pc_file("Thunderegg" "Thunderegg.pc.in" "Thunderegg.pc")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Thunderegg.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
