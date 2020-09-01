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

#ifndef THUNDEREGG_SCHUR_COARSEIFACEINFO_H
#define THUNDEREGG_SCHUR_COARSEIFACEINFO_H
#include <ThunderEgg/PatchInfo.h>
#include <ThunderEgg/Schur/IfaceInfo.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief Represents the interfaces where the neighbor is at a coarser refinement level.
 *
 * There will be two interfaces associated with this object. The interface that lines up with this
 * patch, and interface that lines up with the coarser patch.
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> class CoarseIfaceInfo : public IfaceInfo<D>
{
	public:
	/**
	 * @brief convenience pointer to associated NbrInfo object
	 */
	std::shared_ptr<CoarseNbrInfo<D>> nbr_info;
	/**
	 * @brief The orthant that this patch in relation to the coarser patch's interface.
	 */
	Orthant<D - 1> orth_on_coarse;
	/**
	 * @brief Rank of the coarse interface
	 */
	int coarse_rank;
	/**
	 * @brief The id of the coarser patch's interface
	 */
	int coarse_id;
	/**
	 * @brief The local index of the coarser patch's inteface
	 */
	int coarse_local_index;
	/**
	 * @brief The global index of the coarser patch's interface
	 */
	int coarse_global_index;
	/**
	 * @brief Construct a new CoarseIfaceInfo object
	 *
	 * @param pinfo the cooresponding PatchInfo object
	 * @param s the side that the interface is on
	 */
	CoarseIfaceInfo(std::shared_ptr<PatchInfo<D>> pinfo, Side<D> s)
	{
		nbr_info       = pinfo->getCoarseNbrInfoPtr(s);
		this->id       = pinfo->id * Side<D>::num_sides + s.getIndex();
		orth_on_coarse = nbr_info->orth_on_coarse;
		coarse_id      = nbr_info->id * Side<D>::num_sides + s.opposite().getIndex();
		// fine and coarse interfaces always belong to their patches
		this->rank        = pinfo->rank;
		this->coarse_rank = nbr_info->rank;
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
		ids.push_back(coarse_id);
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
		idx.push_back(coarse_local_index);
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
		idx.push_back(coarse_global_index);
	}
	void getIfaceTypes(std::deque<IfaceType<D>> &types)
	{
		IfaceType<D> fine_type(IfaceType<D>::fine_to_fine, orth_on_coarse);
		IfaceType<D> coarse_type(IfaceType<D>::fine_to_coarse, orth_on_coarse);
		types.push_back(fine_type);
		types.push_back(coarse_type);
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(nbr_info->rank);
		ranks.push_back(nbr_info->rank);
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
		auto it           = rev_map.find(coarse_id);
		if (it != rev_map.end())
			coarse_local_index = it->second;
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index  = rev_map.at(this->local_index);
		coarse_global_index = rev_map.at(coarse_local_index);
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif