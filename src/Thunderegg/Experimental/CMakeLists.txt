if(Zoltan_FOUND)
    list(APPEND Thunderegg_HDRS Thunderegg/Experimental/DomGen.h)
    list(APPEND Thunderegg_SRCS Thunderegg/Experimental/DomGen.cpp)
endif(Zoltan_FOUND)

list(APPEND Thunderegg_HDRS Thunderegg/Experimental/OctNode.h)

list(APPEND Thunderegg_HDRS Thunderegg/Experimental/OctTree.h)

list(APPEND Thunderegg_HDRS Thunderegg/Experimental/PBMatrix.h)
list(APPEND Thunderegg_SRCS Thunderegg/Experimental/PBMatrix.cpp)