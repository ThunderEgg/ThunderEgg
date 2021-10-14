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

#ifndef THUNDEREGG_RUNTIMEERROR_H
#define THUNDEREGG_RUNTIMEERROR_H
/**
 * @file
 *
 * @brief RuntimeError struct
 */

#include <stdexcept>
#include <string>

namespace ThunderEgg
{
/**
 * @brief ThunderEgg runtime exception
 */
struct RuntimeError : std::runtime_error {
	/**
	 * @brief Construct a new RuntimeError object
	 *
	 * @param message the message to print
	 */
	RuntimeError(std::string message) : std::runtime_error(message){};
};
} // namespace ThunderEgg
#endif