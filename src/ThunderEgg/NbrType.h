/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/
#ifndef THUNDEREGG_NBRTYPE_H
#define THUNDEREGG_NBRTYPE_H
/**
 * @file
 *
 * @brief NbrType class
 */
#include <ThunderEgg/tpl/json_fwd.hpp>
#include <ostream>

namespace ThunderEgg {
/**
 * @brief The type of neighbor
 */
enum class NbrType
{
  /**
   * @brief The neighbor is at the same refinement level.
   */
  Normal,
  /**
   * @brief The neighbor is at a coarser refinement level.
   */
  Coarse,
  /**
   * @brief The nighbor is at a finer refinement level.
   */
  Fine
};
/**
 * @brief ostream operator that prints a string representation of NbrType enum.
 */
inline std::ostream&
operator<<(std::ostream& os, const NbrType& type)
{
  switch (type) {
    case NbrType::Coarse:
      os << "NbrType::Coarse";
      break;
    case NbrType::Fine:
      os << "NbrType::Fine";
      break;
    case NbrType::Normal:
      os << "NbrType::Normal";
      break;
  }
  return os;
}

void
to_json(tpl::nlohmann::json& j, const NbrType& o);
void
from_json(const tpl::nlohmann::json& j, NbrType& o);

} // namespace ThunderEgg
#endif