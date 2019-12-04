add_executable(unit_tests 
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
    )
    target_include_directories(unit_tests PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests
    Thunderegg
    )
if(VTK_FOUND)
    target_link_libraries(unit_tests VtkWriter)
endif(VTK_FOUND)
