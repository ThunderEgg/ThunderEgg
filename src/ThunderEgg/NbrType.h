/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
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
#include <ThunderEgg/tpl/json.hpp>
#include <ostream>

namespace ThunderEgg
{
/**
 * @brief The type of neighbor
 */
enum class NbrType {
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
inline std::ostream &operator<<(std::ostream &os, const NbrType &type)
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
NLOHMANN_JSON_SERIALIZE_ENUM(NbrType, {{NbrType::Normal, "NORMAL"}, {NbrType::Coarse, "COARSE"}, {NbrType::Fine, "FINE"}});
} // namespace ThunderEgg
#endif