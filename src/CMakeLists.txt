project(ThunderEgg_Lib)

# determine sources first
list(APPEND ThunderEgg_HDRS ThunderEgg/BiCGStab.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/BiCGStabPatchSolver.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/BiLinearGhostFiller.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/BiLinearGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/BiQuadraticGhostFiller.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/BiQuadraticGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/BufferReader.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/BufferWriter.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/CoarseNbrInfo.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Domain.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Domain.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/DomainGenerator.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/DomainTools.h)

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

list(APPEND ThunderEgg_HDRS ThunderEgg/Operator.h)

if(p4est_FOUND)

  list(APPEND ThunderEgg_HDRS ThunderEgg/P4estDomGen.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/P4estDomGen.cpp)

endif(p4est_FOUND)

list(APPEND ThunderEgg_HDRS ThunderEgg/PatchInfo.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/PatchInfo.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/PatchOperator.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/PatchSolver.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/RuntimeError.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Serializable.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Side.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Timer.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/TriLinearGhostFiller.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/TriLinearGhostFiller.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/TypeDefs.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/ValVector.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/ValVectorGenerator.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Vector.h)
list(APPEND ThunderEgg_SRCS ThunderEgg/Vector.cpp)

list(APPEND ThunderEgg_HDRS ThunderEgg/VectorGenerator.h)

include(ThunderEgg/Experimental/CMakeLists.txt)
include(ThunderEgg/GMG/CMakeLists.txt)
include(ThunderEgg/PETSc/CMakeLists.txt)
include(ThunderEgg/Poisson/CMakeLists.txt)
include(ThunderEgg/Schur/CMakeLists.txt)
include(ThunderEgg/VarPoisson/CMakeLists.txt)

add_library(ThunderEgg ${ThunderEgg_SRCS} ${ThunderEgg_HDRS})
target_compile_features(ThunderEgg PUBLIC cxx_std_11)

# foreach(HDR ${ThunderEgg_HDRS}) get_filename_component(DIR ${HDR} PATH)
# file(COPY ${HDR} DESTINATION include/${DIR}) endforeach()

# set include directories
target_include_directories(
  ThunderEgg PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
list(APPEND ThunderEgg_Includes ${MPI_C_INCLUDE_PATH})
list(APPEND ThunderEgg_Includes ${PETSC_INCLUDES})
list(APPEND ThunderEgg_Includes ${Zoltan_INCLUDES})
if(p4est_FOUND)
  list(APPEND ThunderEgg_Includes ${p4est_INCLUDES})
  list(APPEND ThunderEgg_Includes ${sc_INCLUDES})
endif()
target_include_directories(ThunderEgg PUBLIC ${ThunderEgg_Includes})

# set link list(APPEND ThunderEgg_Libs ${MPI_CXX_LIBRARIES})
list(APPEND ThunderEgg_Libs ${PETSC_LIBRARIES})
if(FFTW_FOUND)
  list(APPEND ThunderEgg_Libs ${FFTW_LIBRARIES})
  list(APPEND ThunderEgg_Includes ${FFTW_INCLUDES})
endif(FFTW_FOUND)
list(APPEND ThunderEgg_Libs ${Zoltan_LIBRARIES})
list(APPEND ThunderEgg_Libs ${BLAS_LIBRARIES})
list(APPEND ThunderEgg_Libs ${LAPACK_LIBRARIES})
if(p4est_FOUND)
  list(APPEND ThunderEgg_Libs ${p4est_LIBRARIES})
  list(APPEND ThunderEgg_Libs ${sc_LIBRARIES})
endif()
target_link_libraries(ThunderEgg PUBLIC ${ThunderEgg_Libs})
target_link_libraries(ThunderEgg PUBLIC ${MPI_C_LIBRARIES})
target_link_libraries(ThunderEgg PUBLIC ${CMAKE_DL_LIBS})

install(
  TARGETS ThunderEgg
  EXPORT ThunderEggConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(
  DIRECTORY ThunderEgg/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*.hpp")

install(EXPORT ThunderEggConfig DESTINATION share/ThunderEgg/cmake)

include(PkgConfigHandler)
configure_pc_file("ThunderEgg" "ThunderEgg.pc.in" "ThunderEgg.pc")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ThunderEgg.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
