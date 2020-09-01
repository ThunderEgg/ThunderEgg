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
template <size_t D> class IfaceInfo
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
	 * @brief the local index in the interface vector.
	 */
	int local_index;
	/**
	 * @brief the global index in the interface vector.
	 */
	int global_index;
	/**
	 * @brief Construct a new IfaceInfo object
	 *
	 * @param rank the rank of the interface
	 * @param id the id of the interface
	 * @param local_index the local index of the interface
	 * @param global_index the global index of the interface
	 */
	IfaceInfo(int rank, int id, int local_index, int global_index)
	: rank(rank), id(id), local_index(local_index), global_index(global_index)
	{
	}
	/**
	 * @brief Destroy the IfaceInfo object
	 */
	virtual ~IfaceInfo() {}
	/**
	 * @brief add to a deque of globally unique ids
	 *
	 * @param ids adds to the deque of ids
	 */
	virtual void getIds(std::deque<int> &ids) = 0;
	/**
	 * @brief add to a deque of local interface indexes
	 *
	 * @param idx adds to the deque of local interface indexes
	 */
	virtual void getLocalIndexes(std::deque<int> &idx) = 0;
	/**
	 * @brief add to a deque of global interface indexes
	 *
	 * @param idx adds to the deque of global interface indexes
	 */
	virtual void getGlobalIndexes(std::deque<int> &idx) = 0;
	/**
	 * @brief add to a deque of IfaceTypes
	 *
	 * @param types adds to the deque of IfaceTypes
	 */
	virtual void getIfaceTypes(std::deque<IfaceType<D>> &types) = 0;
	/**
	 * @brief add to a deque of interface ranks
	 *
	 * @param ranks adds to the deque of interface ranks
	 */
	virtual void getRanks(std::deque<int> &ranks) = 0;
	/**
	 * @brief Set the local indexes in the IfaceInfo objects
	 *
	 * @param rev_map map from id to local_index
	 */
	virtual void setLocalIndexes(const std::map<int, int> &rev_map) = 0;
	/**
	 * @brief Set the global indexes in the IfaceInfo objects
	 *
	 * @param rev_map map form local_index to global_index
	 */
	virtual void setGlobalIndexes(const std::map<int, int> &rev_map) = 0;
};
} // namespace Schur
} // namespace ThunderEgg
#endif