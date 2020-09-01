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

#ifndef THUNDEREGG_SCHUR_NORMALIFACEINFO_H
#define THUNDEREGG_SCHUR_NORMALIFACEINFO_H
#include <ThunderEgg/PatchInfo.h>
#include <ThunderEgg/Schur/IfaceInfo.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief This represents an interface where the neighbor is at the same refinement level
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> class NormalIfaceInfo : public IfaceInfo<D>
{
	public:
	/**
	 * @brief convenience pointer to associated NbrInfo object
	 */
	std::shared_ptr<NormalNbrInfo<D>> nbr_info;
	/**
	 * @brief Construct a new NormalIfaceInfo object all values are set to zero
	 */
	NormalIfaceInfo()
	{
		this->id           = 0;
		this->local_index  = 0;
		this->global_index = 0;
	}
	/**
	 * @brief Construct a new NormalIfaceInfo object
	 *
	 * @param pinfo the associated PatchInfo object
	 * @param s the side of the patch that the interface is on
	 */
	NormalIfaceInfo(std::shared_ptr<PatchInfo<D>> pinfo, Side<D> s)
	{
		nbr_info = pinfo->getNormalNbrInfoPtr(s);
		if (s.isLowerOnAxis()) {
			this->id = pinfo->id * Side<D>::num_sides + s.getIndex();
			// lower axis interface belongs to neighboring rank
			this->rank = nbr_info->rank;
		} else {
			this->id = nbr_info->id * Side<D>::num_sides + s.opposite().getIndex();
			// higher axis interafce belongs to this patch's rank
			this->rank = pinfo->rank;
		}
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
	}
	void getIfaceTypes(std::deque<IfaceType<D>> &types)
	{
		types.push_back(IfaceType<D>::normal);
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(nbr_info->rank);
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index = rev_map.at(this->local_index);
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
