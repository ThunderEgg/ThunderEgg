list(APPEND ThunderEgg_HDRS ThunderEgg/Experimental/OctNode.h)

list(APPEND ThunderEgg_HDRS ThunderEgg/Experimental/OctTree.h)

if(TARGET PETSC::PETSC)
    list(APPEND ThunderEgg_HDRS ThunderEgg/Experimental/PBMatrix.h)
    list(APPEND ThunderEgg_SRCS ThunderEgg/Experimental/PBMatrix.cpp)
endif(TARGET PETSC::PETSC)