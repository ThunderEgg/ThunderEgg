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
template <int D> class NormalIfaceInfo : public IfaceInfo<D>
{
	private:
	/**
	 * @brief Get the rank for the interface on a given side of the patch
	 *
	 * @param pinfo the patch
	 * @param s the side
	 * @return int the rank
	 */
	static int GetRank(std::shared_ptr<const PatchInfo<D>> pinfo, Side<D> s)
	{
		if (s.isLowerOnAxis()) {
			// lower axis interface belongs to neighboring rank
			auto nbr_info = pinfo->getNormalNbrInfo(s);
			return nbr_info.rank;
		} else {
			// higher axis interafce belongs to this patch's rank
			return pinfo->rank;
		}
	}

	/**
	 * @brief Get the id for the interface on a given side of the patch
	 *
	 * @param pinfo the patch
	 * @param s the side
	 * @return int the id
	 */
	static int GetId(std::shared_ptr<const PatchInfo<D>> pinfo, Side<D> s)
	{
		if (s.isLowerOnAxis()) {
			// lower axis interface belongs to neighboring rank
			auto nbr_info = pinfo->getNormalNbrInfo(s);
			return (int) (nbr_info.id * Side<D>::num_sides + s.opposite().getIndex());
		} else {
			// higher axis interafce belongs to this patch's rank
			return (int) (pinfo->id * Side<D>::num_sides + s.getIndex());
		}
	}

	public:
	/**
	 * @brief Construct a new NormalIfaceInfo object
	 *
	 * @param pinfo the associated PatchInfo object
	 * @param s the side of the patch that the interface is on
	 */
	NormalIfaceInfo(std::shared_ptr<const PatchInfo<D>> pinfo, Side<D> s) : IfaceInfo<D>(GetRank(pinfo, s), GetId(pinfo, s)) {}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
