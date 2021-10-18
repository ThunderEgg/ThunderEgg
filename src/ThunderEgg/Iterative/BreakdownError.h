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

#ifndef THUNDEREGG_ITERATIVE_BREAKDOWNERROR_H
#define THUNDEREGG_ITERATIVE_BREAKDOWNERROR_H
/**
 * @file
 *
 * @brief BreakdownError struct
 */

#include <ThunderEgg/RuntimeError.h>

namespace ThunderEgg::Iterative {
/**
 * @brief Breakdown exception for iterative methods
 */
struct BreakdownError : RuntimeError
{
  /**
   * @brief Construct a new RuntimeError object
   *
   * @param message the message to print
   */
  BreakdownError(std::string message)
    : RuntimeError(message){};
};
} // namespace ThunderEgg::Iterative
#endif