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
    Utils.cpp)

target_link_libraries(unit_tests
    Thunderegg
    )
if(VTK_FOUND)
    target_link_libraries(unit_tests VtkWriter)
endif(VTK_FOUND)
