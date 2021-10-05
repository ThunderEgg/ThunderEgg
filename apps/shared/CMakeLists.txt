#determine sources first
if(TARGET P4EST::P4EST)
    list(APPEND UTILS_SRCS GetP8estDomainGenerator.cpp)
    list(APPEND UTILS_SRCS GetP4estDomainGenerator.cpp)
endif()

add_library(UTILS ${UTILS_SRCS})


target_include_directories(UTILS PUBLIC ".")

target_link_libraries(UTILS ThunderEgg)
target_link_libraries(UTILS UTILS_WRITERS)

add_subdirectory(Writers)
