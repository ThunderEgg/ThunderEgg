add_executable(
  unit_tests_mpiany
  test.cpp
  SchurInfo.cpp
  PatchInfo.cpp
  GMG.cpp
  test.cpp
  Side.cpp
  Octant.cpp
  OctTree.cpp
  Domain.cpp
  LocalData.cpp
  DomainTools.cpp
  SimpleGhostFiller.cpp
  BiLinearGhostFiller.cpp
  Utils.cpp
  GMG/CycleFactory.cpp
  GMG/LinearRestrictor.cpp
  utils/DomainReader.cpp
  Timer.cpp)
target_include_directories(unit_tests_mpiany PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests_mpiany Thunderegg)
if(VTK_FOUND)
  target_link_libraries(unit_tests_mpiany VtkWriter)
endif(VTK_FOUND)

add_executable(
  unit_tests_mpi2 test.cpp GMG/InterLevelComm_MPI2.cpp
                  GMG/LinearRestrictor_MPI2.cpp utils/DomainReader.cpp)
target_include_directories(unit_tests_mpi2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests_mpi2 Thunderegg)
if(VTK_FOUND)
  target_link_libraries(unit_tests_mpi2 VtkWriter)
endif(VTK_FOUND)
