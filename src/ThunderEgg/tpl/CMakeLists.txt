# -- add sources

list(APPEND ThunderEgg_tpl_HDRS json.hpp)
list(APPEND ThunderEgg_tpl_HDRS json_fwd.hpp)

# -- install public headers

install(
  FILES ${ThunderEgg_tpl_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg/tpl
  COMPONENT   ThunderEgg_Development
)
