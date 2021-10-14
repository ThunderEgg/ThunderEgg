/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

#ifndef THUNDEREGG_GHOSTFILLINGTYPE_H
#define THUNDEREGG_GHOSTFILLINGTYPE_H
/**
 * @file
 *
 * @brief GhostFillingType enum
 */

namespace ThunderEgg
{
/**
 * @brief type of ghost filling.
 */
enum class GhostFillingType {
	/**
	 * @brief Fill only faces
	 */
	Faces,
	/**
	 * @brief Fill faces and edges
	 */
	Edges,
	/**
	 * @brief Fill faces, edges, and corners. (or faces and corners in the 2d case)
	 */
	Corners
};
} // namespace ThunderEgg

#endif