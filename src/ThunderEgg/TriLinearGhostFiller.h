/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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
/**
 * @file
 *
 * @brief TriLinearGhostFiller class
 */

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
	void fillGhostCellsForNbrPatch(const PatchInfo<3> &              pinfo,
	                               const PatchView<const double, 3> &local_view,
	                               const PatchView<const double, 3> &nbr_view,
	                               Side<3>                           side,
	                               NbrType                           nbr_type,
	                               Orthant<2>                        orthant_on_coarse) const override;

	void fillGhostCellsForEdgeNbrPatch(const PatchInfo<3> &              pinfo,
	                                   const PatchView<const double, 3> &local_view,
	                                   const PatchView<const double, 3> &nbr_view,
	                                   Edge                              edge,
	                                   NbrType                           nbr_type,
	                                   Orthant<1>                        orthant_on_coarse) const override;

	void fillGhostCellsForCornerNbrPatch(const PatchInfo<3> &              pinfo,
	                                     const PatchView<const double, 3> &local_view,
	                                     const PatchView<const double, 3> &nbr_view,
	                                     Corner<3>                         corner,
	                                     NbrType                           nbr_type) const override;

	void fillGhostCellsForLocalPatch(const PatchInfo<3> &pinfo, const PatchView<const double, 3> &view) const override;
	/**
	 * @brief Construct a new TriLinearGhostFiller object
	 *
	 * Currently, this only supports an even number of cells on each axis of the patch
	 *
	 * @param domain the domain on which ghosts will be filled
	 * @param fill_type the ghost filling type to perform
	 */
	TriLinearGhostFiller(const Domain<3> &domain, GhostFillingType fill_type);
	/**
	 * @brief Clone this TriLienarGhostFiller
	 *
	 * @return TriLinearGhostFiller* a newly allocated copy
	 */
	TriLinearGhostFiller *clone() const override;
};
} // namespace ThunderEgg
#endif
