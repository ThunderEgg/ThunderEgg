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
/**
 * @file
 *
 * @brief CoarseIfaceInfo class
 */
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
template <int D> class CoarseIfaceInfo : public IfaceInfo<D>
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
	 * @brief The local column index of the coarser patch's inteface
	 */
	int coarse_col_local_index = -1;
	/**
	 * @brief The global index of the coarser patch's interface
	 */
	int coarse_global_index = -1;
	/**
	 * @brief Construct a new CoarseIfaceInfo object
	 *
	 * indexes will be set to -1
	 *
	 * @param pinfo the cooresponding PatchInfo object
	 * @param s the side that the interface is on
	 */
	CoarseIfaceInfo(const PatchInfo<D> &pinfo, Side<D> s) : IfaceInfo<D>(pinfo.rank, GetId(pinfo, s))
	{
		// fine and coarse interfaces always belong to their patches
		auto nbr_info  = pinfo.getCoarseNbrInfo(s);
		orth_on_coarse = nbr_info.orth_on_coarse;
		coarse_rank    = nbr_info.rank;
		coarse_id      = nbr_info.id * Side<D>::number_of + s.opposite().getIndex();
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif