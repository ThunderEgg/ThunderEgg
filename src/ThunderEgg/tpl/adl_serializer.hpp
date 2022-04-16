#pragma once

#include <type_traits>
#include <utility>

#include <ThunderEgg/tpl/detail/conversions/from_json.hpp>
#include <ThunderEgg/tpl/detail/conversions/to_json.hpp>
#include <ThunderEgg/tpl/detail/meta/identity_tag.hpp>
#include <ThunderEgg/tpl/detail/meta/type_traits.hpp>

namespace ThunderEgg::tpl::nlohmann {

/// @sa https://json.nlohmann.me/api/adl_serializer/
template<typename ValueType, typename>
struct adl_serializer
{
  /// @brief convert a JSON value to any value type
  /// @sa https://json.nlohmann.me/api/adl_serializer/from_json/
  template<typename BasicJsonType, typename TargetType = ValueType>
  static auto from_json(BasicJsonType&& j, TargetType& val) noexcept(
    noexcept(::ThunderEgg::tpl::nlohmann::from_json(std::forward<BasicJsonType>(j), val)))
    -> decltype(::ThunderEgg::tpl::nlohmann::from_json(std::forward<BasicJsonType>(j), val), void())
  {
    ::ThunderEgg::tpl::nlohmann::from_json(std::forward<BasicJsonType>(j), val);
  }

  /// @brief convert a JSON value to any value type
  /// @sa https://json.nlohmann.me/api/adl_serializer/from_json/
  template<typename BasicJsonType, typename TargetType = ValueType>
  static auto from_json(BasicJsonType&& j) noexcept(
    noexcept(::ThunderEgg::tpl::nlohmann::from_json(std::forward<BasicJsonType>(j),
                                                    detail::identity_tag<TargetType>{})))
    -> decltype(::ThunderEgg::tpl::nlohmann::from_json(std::forward<BasicJsonType>(j),
                                                       detail::identity_tag<TargetType>{}))
  {
    return ::ThunderEgg::tpl::nlohmann::from_json(std::forward<BasicJsonType>(j),
                                                  detail::identity_tag<TargetType>{});
  }

  /// @brief convert any value type to a JSON value
  /// @sa https://json.nlohmann.me/api/adl_serializer/to_json/
  template<typename BasicJsonType, typename TargetType = ValueType>
  static auto to_json(BasicJsonType& j, TargetType&& val) noexcept(
    noexcept(::ThunderEgg::tpl::nlohmann::to_json(j, std::forward<TargetType>(val))))
    -> decltype(::ThunderEgg::tpl::nlohmann::to_json(j, std::forward<TargetType>(val)), void())
  {
    ::ThunderEgg::tpl::nlohmann::to_json(j, std::forward<TargetType>(val));
  }
};
} // namespace ThunderEgg::tpl::nlohmann
