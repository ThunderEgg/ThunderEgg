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

#ifndef THUNDEREGG_SCHUR_IFACEINFO_H
#define THUNDEREGG_SCHUR_IFACEINFO_H
#include <ThunderEgg/Schur/IfaceType.h>
#include <deque>
#include <map>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief The IfaceInfo class represents the information for an interface on a given side of the
 * patch.
 *
 * The information contained  will be the globally unique ID and the local and global index(es) in
 * the interface vector.
 *
 * @tparam D the number of Cartesian dimensions in the patches.
 */
template <int D> class IfaceInfo
{
	public:
	/**
	 * @brief The rank that the interface resides on.
	 */
	int rank;
	/**
	 * @brief The globally unique ID of the interface.
	 */
	int id;
	/**
	 * @brief the local index in the interface for the local patch interface vector
	 */
	int patch_local_index = -1;
	/**
	 * @brief the local index in the interface for the local column interface vector
	 */
	int col_local_index = -1;
	/**
	 * @brief the local index in the interface for the local row interface vector
	 */
	int row_local_index = -1;
	/**
	 * @brief the global index in the interface vector.
	 */
	int global_index = -1;
	/**
	 * @brief Construct a new IfaceInfo object
	 *
	 * All indexes will be set to -1
	 *
	 * @param rank the rank of the interface
	 * @param id the id of the interface
	 */
	IfaceInfo(int rank, int id) : rank(rank), id(id) {}
	/**
	 * @brief Destroy the IfaceInfo object
	 */
	virtual ~IfaceInfo() {}
};
} // namespace Schur
} // namespace ThunderEgg
#endif