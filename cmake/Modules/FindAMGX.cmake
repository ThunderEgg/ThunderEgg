# - Try to find AMGX
#

find_path (AMGX_DIR NAMES include/amgx_c.h HINTS ENV AMGX_DIR)

IF(AMGX_DIR)
  SET(AMGX_INCLUDES ${AMGX_DIR})
  find_path (AMGX_INCLUDE_DIR NAMES amgx_c.h HINTS "${AMGX_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND AMGX_INCLUDES ${AMGX_INCLUDE_DIR})
  FILE(GLOB AMGX_LIBRARIES RELATIVE "${AMGX_DIR}/lib" "${AMGX_DIR}/lib/libamgxsh.so")
ENDIF(AMGX_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMGX REQUIRED_VARS AMGX_LIBRARIES AMGX_INCLUDE_DIR)
