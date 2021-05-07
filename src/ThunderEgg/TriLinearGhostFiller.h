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

#ifndef THUNDEREGG_TRILINEARGHOSTFILLER_H
#define THUNDEREGG_TRILINEARGHOSTFILLER_H
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/MPIGhostFiller.h>
namespace ThunderEgg
{
/**
 * @brief Performs trilinear interpolation on coarse-fine boundaries of patches
 *
 * It only uses the coarse cell, and the four cooresponding fine cells to interpolate on the
 * coarse-fine boundary.
 */
class TriLinearGhostFiller : public MPIGhostFiller<3>
{
	public:
	void fillGhostCellsForNbrPatch(const PatchInfo<3> &             pinfo,
	                               const std::vector<LocalData<3>> &local_datas,
	                               const std::vector<LocalData<3>> &nbr_datas,
	                               Side<3>                          side,
	                               NbrType                          nbr_type,
	                               Orthant<2>                       orthant_on_coarse) const override;

	void fillGhostCellsForEdgeNbrPatch(const PatchInfo<3> &             pinfo,
	                                   const std::vector<LocalData<3>> &local_datas,
	                                   const std::vector<LocalData<3>> &nbr_datas,
	                                   Edge<3>                          edge,
	                                   NbrType                          nbr_type,
	                                   Orthant<1>                       orthant_on_coarse) const override
	{
	}

	void fillGhostCellsForCornerNbrPatch(const PatchInfo<3> &             pinfo,
	                                     const std::vector<LocalData<3>> &local_datas,
	                                     const std::vector<LocalData<3>> &nbr_datas,
	                                     Corner<3>                        corner,
	                                     NbrType                          nbr_type) const override
	{
	}

	void fillGhostCellsForLocalPatch(const PatchInfo<3> &pinfo, const std::vector<LocalData<3>> &local_datas) const override;
	/**
	 * @brief Construct a new TriLinearGhostFiller object
	 *
	 * Currently, this only supports an even number of cells on each axis of the patch
	 *
	 * @param domain the domain on which ghosts will be filled
	 */
	explicit TriLinearGhostFiller(std::shared_ptr<const Domain<3>> domain);
};
} // namespace ThunderEgg
#endif
