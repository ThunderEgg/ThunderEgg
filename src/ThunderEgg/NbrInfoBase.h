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
#ifndef THUNDEREGG_NBRINFOBASE_H
#define THUNDEREGG_NBRINFOBASE_H
/**
 * @file
 *
 * @brief NbrInfoBase class
 */
#include <ThunderEgg/NbrType.h>
#include <ThunderEgg/Serializable.h>
#include <deque>
#include <map>
#include <memory>

namespace ThunderEgg
{
/**
 * @brief Represents information about a patch's neighbor.
 *
 * Includes information like neighbor id and and indexes.
 */
class NbrInfoBase : public Serializable
{
	public:
	/**
	 * @brief Destroy the NbrInfo object
	 */
	virtual ~NbrInfoBase() = default;
	/**
	 * @brief Get the NbrType
	 */
	virtual NbrType getNbrType() const = 0;
	/**
	 * @brief Add to a deque of neighbor ids
	 */
	virtual void getNbrIds(std::deque<int> &nbr_ids) const = 0;
	/**
	 * @brief Add to a deque of neighbor ranks
	 */
	virtual void getNbrRanks(std::deque<int> &nbr_ranks) const = 0;
	/**
	 * @brief Set the local indexes in the NbrInfo objects
	 *
	 * @param rev_map map from id to local_index
	 */
	virtual void setGlobalIndexes(const std::map<int, int> &rev_map) = 0;
	/**
	 * @brief Set the global indexes in the NbrInfo objects
	 *
	 * @param rev_map map from local_index to global_index
	 */
	virtual void setLocalIndexes(const std::map<int, int> &rev_map) = 0;
	/**
	 * @brief get a clone of this object (equivalent to copy constructor)
	 *
	 * @return std::unique_ptr<NbrInfo> the cloned object
	 */
	virtual std::unique_ptr<NbrInfoBase> clone() const = 0;
};
} // namespace ThunderEgg
#endif