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

#ifndef THUNDEREGG_SCHUR_FINEIFACEINFO_H
#define THUNDEREGG_SCHUR_FINEIFACEINFO_H
#include <ThunderEgg/PatchInfo.h>
#include <ThunderEgg/Schur/IfaceInfo.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief Represents the interfaces where the neighbors are at a finer refinement level.
 *
 * There will be 2^(D-1)+1 interfaces associated with this object. The interface that lines up with
 * this patch, and interface that lines up with the finer patches.
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> class FineIfaceInfo : public IfaceInfo<D>
{
	public:
	/**
	 * @brief convenience pointer to associated NbrInfo object
	 */
	std::shared_ptr<FineNbrInfo<D>> nbr_info;
	/**
	 * @brief the ranks of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_ranks;
	/**
	 * @brief the ids of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_ids;
	/**
	 * @brief the local indexes of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_local_indexes;
	/**
	 * @brief the global indexes of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_global_indexes;
	/**
	 * @brief Construct a new FineIfaceInfo object
	 *
	 * @param pinfo the associated PatchInfo object
	 * @param s the side of the patch that the interface is on
	 */
	FineIfaceInfo(std::shared_ptr<PatchInfo<D>> pinfo, Side<D> s)
	{
		nbr_info   = pinfo->getFineNbrInfoPtr(s);
		this->id   = pinfo->id * Side<D>::num_sides + s.getIndex();
		this->rank = pinfo->rank;
		for (size_t i = 0; i < fine_ids.size(); i++) {
			fine_ids[i]   = nbr_info->ids[i] * Side<D>::num_sides + s.opposite().getIndex();
			fine_ranks[i] = nbr_info->ranks[i];
		}
	}
	void getIdxAndTypes(std::deque<int> &idx, std::deque<IfaceType<D>> &types)
	{
		idx.push_back(this->local_index);
		types.push_back(IfaceType<D>::coarse_to_coarse);
		for (size_t i = 0; i < fine_local_indexes.size(); i++) {
			idx.push_back(fine_local_indexes[i]);
			IfaceType<D> type(IfaceType<D>::coarse_to_fine, i);
			types.push_back(type);
		}
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			ids.push_back(fine_ids[i]);
		}
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			idx.push_back(fine_local_indexes[i]);
		}
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			idx.push_back(fine_global_indexes[i]);
		}
	}
	void getIfaceTypes(std::deque<IfaceType<D>> &types)
	{
		types.push_back(IfaceType<D>::coarse_to_coarse);
		for (Orthant<D - 1> o : Orthant<D - 1>::getValues()) {
			IfaceType<D> type(IfaceType<D>::coarse_to_fine, o);
			types.push_back(type);
		}
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(-1);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			ranks.push_back(nbr_info->ranks[i]);
		}
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			auto it = rev_map.find(this->fine_ids[i]);
			if (it != rev_map.end())
				fine_local_indexes[i] = it->second;
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index = rev_map.at(this->local_index);
		for (size_t i = 0; i < fine_local_indexes.size(); i++) {
			fine_global_indexes[i] = rev_map.at(fine_local_indexes[i]);
		}
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif