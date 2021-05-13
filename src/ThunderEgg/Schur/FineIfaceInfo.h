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
template <int D> class FineIfaceInfo : public IfaceInfo<D>
{
	private:
	/**
	 * @brief Get the id for the interface on a given side of the patch
	 *
	 * @param pinfo the patch
	 * @param s the side
	 * @return int the id
	 */
	static int GetId(const PatchInfo<D> &pinfo, Side<D> s)
	{
		return (int) (pinfo.id * Side<D>::number_of + s.getIndex());
	}

	public:
	/**
	 * @brief the ranks of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_ranks;
	/**
	 * @brief the ids of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_ids;
	/**
	 * @brief the local column indexes of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_col_local_indexes;
	/**
	 * @brief the global indexes of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_global_indexes;
	/**
	 * @brief Construct a new FineIfaceInfo object
	 *
	 * indexes will be set to -1
	 *
	 * @param pinfo the associated PatchInfo object
	 * @param s the side of the patch that the interface is on
	 */
	FineIfaceInfo(const PatchInfo<D> &pinfo, Side<D> s) : IfaceInfo<D>(pinfo.rank, GetId(pinfo, s))
	{
		auto nbr_info = pinfo.getFineNbrInfo(s);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			fine_ids[i]   = nbr_info.ids[i] * Side<D>::number_of + s.opposite().getIndex();
			fine_ranks[i] = nbr_info.ranks[i];
		}
		fine_col_local_indexes.fill(-1);
		fine_global_indexes.fill(-1);
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif