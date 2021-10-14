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

#include <ThunderEgg/BiQuadraticGhostFiller.h>
#include <ThunderEgg/RuntimeError.h>
namespace ThunderEgg
{
BiQuadraticGhostFiller::BiQuadraticGhostFiller(const Domain<2> &domain, GhostFillingType fill_type) : MPIGhostFiller<2>(domain, fill_type) {}
BiQuadraticGhostFiller *BiQuadraticGhostFiller::clone() const
{
	return new BiQuadraticGhostFiller(*this);
}
namespace
{
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a coarser neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForFineNbrLower and FillGhostForFineNbrUpper
 *
 * @param view the patch data
 * @param side the side that the neighbor patch is on
 */
void FillGhostForLocalWithCoarseNbr(const PatchView<const double, 2> &view, Side<2> side)
{
	View<const double, 2> inner_slice = view.getSliceOn(side, {1});
	View<const double, 2> slice       = view.getSliceOn(side, {0});
	View<double, 2>       ghost       = view.getGhostSliceOn(side, {0});
	int                   n           = ghost.getEnd()[0] + 1;
	for (int c = slice.getStart()[1]; c <= slice.getEnd()[1]; c++) {
		for (int idx = 0; idx < n; idx++) {
			ghost(idx, c) += 2 * slice(idx, c) / 3 - inner_slice(idx, c) / 5;
		}
	}
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a finer neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForCoarseNbrLower and FillGhostForCoarseNbrUpper
 *
 * @param view the patch data
 * @param side the side that the neighbor patch is on
 */
void FillGhostForLocalWithFineNbr(const PatchView<const double, 2> &view, Side<2> side)
{
	View<const double, 2> slice = view.getSliceOn(side, {0});
	View<double, 2>       ghost = view.getGhostSliceOn(side, {0});
	int                   n     = ghost.getEnd()[0] + 1;
	for (int c = slice.getStart()[1]; c <= slice.getEnd()[1]; c++) {
		ghost(0, c) += -slice(0, c) / 10 + slice(1, c) / 15 - slice(2, c) / 30;
		for (int idx = 1; idx < n - 1; idx++) {
			ghost(idx, c) += -slice(idx - 1, c) / 30 - slice(idx + 1, c) / 30;
		}
		ghost(n - 1, c) += -slice(n - 1, c) / 10 + slice(n - 2, c) / 15 - slice(n - 3, c) / 30;
	}
}
/**
 * @brief This is just a simple copy of values
 *
 * @param local_view the local pach
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 */
void FillGhostForNormalNbr(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Side<2> side)
{
	View<const double, 2> local_slice = local_view.getSliceOn(side, {0});
	View<double, 2>       nbr_ghosts  = nbr_view.getGhostSliceOn(side.opposite(), {0});
	Loop::OverInteriorIndexes<2>(nbr_ghosts, [&](const std::array<int, 2> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
}
/**
 * @brief Fill the ghost values for a coarse neighbor when this patch is on the lower part of the coarser neighbors side
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be needed.
 *
 * Those values are added in FillGhostForLocalWithFineNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the coarser neighbor face that this patch lies on
 */
void FillGhostForCoarseNbrLower(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Side<2> side)
{
	View<const double, 2> slice       = local_view.getSliceOn(side, {0});
	View<const double, 2> inner_slice = local_view.getSliceOn(side, {1});
	View<double, 2>       ghost       = nbr_view.getGhostSliceOn(side.opposite(), {0});
	for (int c = slice.getStart()[1]; c <= slice.getEnd()[1]; c++) {
		for (int idx = slice.getStart()[0]; idx <= slice.getEnd()[0]; idx++) {
			ghost(idx / 2, c) += slice(idx, c) / 3 + inner_slice(idx, c) / 5;
		}
	}
}
/**
 * @brief Fill the ghost values for a coarse neighbor when this patch is on the upper part of the coarser neighbors side
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be needed.
 *
 * Those values are added in FillGhostForLocalWithFineNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the coarser neighbor face that this patch lies on
 */
void FillGhostForCoarseNbrUpper(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Side<2> side)
{
	View<const double, 2> slice       = local_view.getSliceOn(side, {0});
	View<const double, 2> inner_slice = local_view.getSliceOn(side, {1});
	View<double, 2>       ghost       = nbr_view.getGhostSliceOn(side.opposite(), {0});
	int                   n           = ghost.getEnd()[0] + 1;
	for (int c = slice.getStart()[1]; c <= slice.getEnd()[1]; c++) {
		for (int idx = slice.getStart()[0]; idx <= slice.getEnd()[0]; idx++) {
			ghost((idx + n) / 2, c) += slice(idx, c) / 3 + inner_slice(idx, c) / 5;
		}
	}
}
/**
 * @brief Fill the ghost values for a fine neighbor for when the neighbor is on the lower part of this patches side
 *
 * This is only part of the value. Values from the interior of the fine neighbor will also be needed.
 *
 * Those values are added in FillGhostForLocalWithCoarseNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the this patches face that the finer neighbor patch lies on
 */
void FillGhostForFineNbrLower(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Side<2> side)
{
	View<const double, 2> slice = local_view.getSliceOn(side, {0});
	View<double, 2>       ghost = nbr_view.getGhostSliceOn(side.opposite(), {0});
	for (int c = slice.getStart()[1]; c <= slice.getEnd()[1]; c++) {
		ghost(0, c) += 3 * slice(0, c) / 4 - 3 * slice(1, c) / 10 + slice(2, c) / 12;
		ghost(1, c) += 7 * slice(0, c) / 20 + 7 * slice(1, c) / 30 - slice(2, c) / 20;
		for (int idx = slice.getStart()[0] + 2; idx <= slice.getEnd()[0]; idx++) {
			if (idx % 2 == 0) {
				ghost(idx, c) += slice(idx / 2 - 1, c) / 12 + slice(idx / 2, c) / 2 - slice(idx / 2 + 1, c) / 20;
			} else {
				ghost(idx, c) += -slice(idx / 2 - 1, c) / 20 + slice(idx / 2, c) / 2 + slice(idx / 2 + 1, c) / 12;
			}
		}
	}
}
/**
 * @brief Fill the ghost values for a fine neighbor for when the neighbor is on the upper part of this patches side
 *
 * This is only part of the value. Values from the interior of the fine neighbor will also be needed.
 *
 * Those values are added in FillGhostForLocalWithCoarseNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the this patches face that the finer neighbor patch lies on
 */
void FillGhostForFineNbrUpper(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Side<2> side)
{
	View<const double, 2> slice = local_view.getSliceOn(side, {0});
	View<double, 2>       ghost = nbr_view.getGhostSliceOn(side.opposite(), {0});
	int                   n     = ghost.getEnd()[0] + 1;
	for (int c = slice.getStart()[1]; c <= slice.getEnd()[1]; c++) {
		for (int idx = 0; idx < n - 2; idx++) {
			if ((idx + n) % 2 == 0) {
				ghost(idx, c) += slice((idx + n) / 2 - 1, c) / 12 + slice((idx + n) / 2, c) / 2 - slice((idx + n) / 2 + 1, c) / 20;
			} else {
				ghost(idx, c) += -slice((idx + n) / 2 - 1, c) / 20 + slice((idx + n) / 2, c) / 2 + slice((idx + n) / 2 + 1, c) / 12;
			}
		}
		ghost(n - 2, c) += 7 * slice(n - 1, c) / 20 + 7 * slice(n - 2, c) / 30 - slice(n - 3, c) / 20;
		ghost(n - 1, c) += 3 * slice(n - 1, c) / 4 - 3 * slice(n - 2, c) / 10 + slice(n - 3, c) / 12;
	}
}

/**
 * @brief This is just a simple copy of values
 *
 * @param local_view the local pach
 * @param nbr_view the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForNormalCornerNbr(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Corner<2> corner)
{
	View<const double, 1> local_slice = local_view.getSliceOn(corner, {0, 0});
	View<double, 1>       nbr_ghosts  = nbr_view.getGhostSliceOn(corner.opposite(), {0, 0});
	for (int c = local_slice.getStart()[0]; c <= local_slice.getEnd()[0]; c++) {
		nbr_ghosts(c) = local_slice(c);
	}
}
/**
 * @brief Fill the ghost values for a coarse neighbor
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be needed.
 *
 * Those values are added in FillGhostForLocalWithCornerFineNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForCoarseCornerNbr(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Corner<2> corner)
{
	View<const double, 1> slice       = local_view.getSliceOn(corner, {0, 0});
	View<const double, 1> inner_slice = local_view.getSliceOn(corner, {1, 1});
	View<double, 1>       ghost       = nbr_view.getGhostSliceOn(corner.opposite(), {0, 0});
	for (int c = slice.getStart()[0]; c <= slice.getEnd()[0]; c++) {
		ghost(c) += 2 * slice(c) / 3 + 2 * inner_slice(c) / 5;
	}
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a finer neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForCoarseCornerNbr
 *
 * @param view the patch data
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForLocalWithFineCornerNbr(const PatchView<const double, 2> &view, Corner<2> corner)
{
	View<const double, 1> slice = view.getSliceOn(corner, {0, 0});
	View<double, 1>       ghost = view.getGhostSliceOn(corner, {0, 0});
	for (int c = slice.getStart()[0]; c <= slice.getEnd()[0]; c++) {
		ghost(c) += -slice(c) / 15;
	}
}

/**
 * @brief Fill the ghost values for a fine neighbor
 *
 * This is only part of the value. Values from the interior of the fine neighbor will also be needed.
 *
 * Those values are added in FillGhostForLocalWithCornerCoarseNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForFineCornerNbr(const PatchView<const double, 2> &local_view, const PatchView<const double, 2> &nbr_view, Corner<2> corner)
{
	View<const double, 1> slice = local_view.getSliceOn(corner, {0, 0});
	View<double, 1>       ghost = nbr_view.getGhostSliceOn(corner.opposite(), {0, 0});
	for (int c = slice.getStart()[0]; c <= slice.getEnd()[0]; c++) {
		ghost(c) += 8 * slice(c) / 15;
	}
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a coarser neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForFineCornerNbr
 *
 * @param view the patch data
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForLocalWithCoarseCornerNbr(const PatchView<const double, 2> &view, Corner<2> corner)
{
	View<const double, 1> inner_slice = view.getSliceOn(corner, {1, 1});
	View<const double, 1> slice       = view.getSliceOn(corner, {0, 0});
	View<double, 1>       ghost       = view.getGhostSliceOn(corner, {0, 0});
	for (int c = slice.getStart()[0]; c <= slice.getEnd()[0]; c++) {
		ghost(c) += 2 * slice(c) / 3 - inner_slice(c) / 5;
	}
}
/**
 * @brief Add in local information to the ghost cells on the sides
 *
 * @param pinfo the patch
 * @param view the patch data
 */
void FillLocalGhostCellsOnSides(const PatchInfo<2> &pinfo, const PatchView<const double, 2> &view)
{
	for (Side<2> side : Side<2>::getValues()) {
		if (pinfo.hasNbr(side)) {
			switch (pinfo.getNbrType(side)) {
				case NbrType::Normal:
					// nothing need to be done
					break;
				case NbrType::Coarse:
					FillGhostForLocalWithCoarseNbr(view, side);
					break;
				case NbrType::Fine:
					FillGhostForLocalWithFineNbr(view, side);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
/**
 * @brief Add in local information to the ghost cells on the corners
 *
 * @param pinfo the patch
 * @param view the patch data
 */
void FillLocalGhostCellsOnCorners(const PatchInfo<2> &pinfo, const PatchView<const double, 2> &view)
{
	for (Corner<2> corner : Corner<2>::getValues()) {
		if (pinfo.hasNbr(corner)) {
			switch (pinfo.getNbrType(corner)) {
				case NbrType::Normal:
					// nothing need to be done
					break;
				case NbrType::Coarse:
					FillGhostForLocalWithCoarseCornerNbr(view, corner);
					break;
				case NbrType::Fine:
					FillGhostForLocalWithFineCornerNbr(view, corner);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
} // namespace

void BiQuadraticGhostFiller::fillGhostCellsForNbrPatch(const PatchInfo<2>               &pinfo,
                                                       const PatchView<const double, 2> &local_view,
                                                       const PatchView<const double, 2> &nbr_view,
                                                       Side<2>                           side,
                                                       NbrType                           nbr_type,
                                                       Orthant<1>                        orthant) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostForNormalNbr(local_view, nbr_view, side);
			break;
		case NbrType::Coarse:
			if (orthant == Orthant<1>::lower()) {
				FillGhostForCoarseNbrLower(local_view, nbr_view, side);
			} else {
				FillGhostForCoarseNbrUpper(local_view, nbr_view, side);
			}
			break;
		case NbrType::Fine:
			if (orthant == Orthant<1>::lower()) {
				FillGhostForFineNbrLower(local_view, nbr_view, side);
			} else {
				FillGhostForFineNbrUpper(local_view, nbr_view, side);
			}
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void BiQuadraticGhostFiller::fillGhostCellsForEdgeNbrPatch(const PatchInfo<2>               &pinfo,
                                                           const PatchView<const double, 2> &local_view,
                                                           const PatchView<const double, 2> &nbr_view,
                                                           Edge                              edge,
                                                           NbrType                           nbr_type,
                                                           Orthant<1>                        orthant_on_coarse) const
{
	// no edges for 2d
}

void BiQuadraticGhostFiller::fillGhostCellsForCornerNbrPatch(const PatchInfo<2>               &pinfo,
                                                             const PatchView<const double, 2> &local_view,
                                                             const PatchView<const double, 2> &nbr_view,
                                                             Corner<2>                         corner,
                                                             NbrType                           nbr_type) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostForNormalCornerNbr(local_view, nbr_view, corner);
			break;
		case NbrType::Coarse:
			FillGhostForCoarseCornerNbr(local_view, nbr_view, corner);
			break;
		case NbrType::Fine:
			FillGhostForFineCornerNbr(local_view, nbr_view, corner);
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void BiQuadraticGhostFiller::fillGhostCellsForLocalPatch(const PatchInfo<2> &pinfo, const PatchView<const double, 2> &view) const
{
	switch (this->getFillType()) {
		case GhostFillingType::Corners: // Fill corners and faces
			FillLocalGhostCellsOnCorners(pinfo, view);
			[[fallthrough]];
		case GhostFillingType::Faces:
			FillLocalGhostCellsOnSides(pinfo, view);
			break;
		default:
			throw RuntimeError("Unsupported GhostFillingType");
	}
}
} // namespace ThunderEgg
