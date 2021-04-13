#determine sources first
list(APPEND UTILS_SRCS Init.cpp Init.h)
if(TARGET P4EST::P4EST)
    list(APPEND UTILS_SRCS TreeToP4est.cpp TreeToP4est.h)
endif()
if(TARGET P4EST::P8EST)
    list(APPEND UTILS_SRCS TreeToP8est.cpp TreeToP8est.h)
endif()

add_library(UTILS ${UTILS_SRCS})


target_include_directories(UTILS PUBLIC ".")

target_link_libraries(UTILS ThunderEgg)
target_link_libraries(UTILS UTILS_WRITERS)

add_subdirectory(Writers)
