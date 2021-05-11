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

#ifndef THUNDEREGG_GHOSTFILLINGTYPE_H
#define THUNDEREGG_GHOSTFILLINGTYPE_H

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