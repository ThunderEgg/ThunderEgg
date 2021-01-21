# - Try to find sc
#

find_path (sc_DIR NAMES include/sc.h HINTS p4est_DIR sc_DIR ENV p4est_DIR ENV sc_DIR CPATH)

  SET(sc_INCLUDES ${sc_DIR})
  find_path (sc_INCLUDE_DIR NAMES sc.h HINTS "${sc_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND sc_INCLUDES ${sc_INCLUDE_DIR})
  find_library(sc_LIBRARIES NAMES sc PATHS "${sc_DIR}/lib" ${sc_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(sc REQUIRED_VARS sc_LIBRARIES sc_INCLUDE_DIR)
