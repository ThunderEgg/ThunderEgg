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

#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/TriLinearGhostFiller.h>
namespace ThunderEgg
{
namespace
{
/**
 * @brief Get the Offset on a given orthant
 *
 * @param end coordnate of non-ghost celsl
 * @param orth the orthant
 * @return std::array<int,2> the offsets
 */
std::array<int, 2> getOffset(const std::array<int, 4> end, Side<3> s, Orthant<2> orth)
{
	std::array<int, 2> offset = {0, 0};
	for (size_t i = 0; i < s.getAxisIndex(); i++) {
		if (orth.isHigherOnAxis(i)) {
			offset[i] = end[i] + 1;
		}
	}
	for (size_t i = s.getAxisIndex() + 1; i < 3; i++) {
		if (orth.isHigherOnAxis(i - 1)) {
			offset[i - 1] = end[i] + 1;
		}
	}
	return offset;
}

///////
// SIDES
///////

/**
 * @brief Simple copy of values
 *
 * @param local_view the neighbor data
 * @param nbr_view  the local data
 * @param side the side that the neighbor is on
 */
void FillGhostCellsForNormalNbr(const PatchView<const double, 3> &local_view, const PatchView<const double, 3> &nbr_view, Side<3> side)
{
	View<const double, 3> local_slice = local_view.getSliceOn(side, {0});
	View<double, 3>       nbr_ghosts  = nbr_view.getGhostSliceOn(side.opposite(), {0});
	loop_over_interior_indexes<3>(nbr_ghosts, [&](const std::array<int, 3> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
}
/**
 * @brief Fill ghost cells for coarser neighbor
 *
 * This is only part of the interpolation, values are also needed from the neighbor patch.
 * These values are added in FillGhostCellsForLocalWithFineNbr
 *
 * @param local_view local patch data
 * @param nbr_view neighbor patch data
 * @param side the side of the patch that the neighbor patch is on
 * @param orthant the orthant of the neighbors side that this patch lies on
 */
void FillGhostCellsForCoarseNbr(const PatchView<const double, 3> &local_view,
                                const PatchView<const double, 3> &nbr_view,
                                Side<3>                           side,
                                Orthant<2>                        orthant)
{
	auto [offset_i, offset_j]         = getOffset(local_view.getEnd(), side, orthant);
	View<const double, 3> local_slice = local_view.getSliceOn(side, {0});
	View<double, 3>       nbr_ghosts  = nbr_view.getGhostSliceOn(side.opposite(), {0});
	for (int c = local_slice.getStart()[2]; c <= local_slice.getEnd()[2]; c++) {
		for (int j = local_slice.getStart()[1]; j <= local_slice.getEnd()[1]; j++) {
			for (int i = local_slice.getStart()[0]; i <= local_slice.getEnd()[0]; i++) {
				nbr_ghosts((i + offset_i) / 2, (j + offset_j) / 2, c) += 1.0 / 3.0 * local_slice(i, j, c);
			}
		}
	}
}
/**
 * @brief Add in the values need for this patch's ghost cells when there is a fine neighbor
 *
 * These values compliment FIllGhostCellsForCoarseNbr
 *
 * @param view the patch view
 * @param side the side of the patch that the neighbor is on
 */
void FillGhostCellsForLocalWithFineNbr(const PatchView<const double, 3> &view, Side<3> side)
{
	View<const double, 3> local_slice  = view.getSliceOn(side, {0});
	View<double, 3>       local_ghosts = view.getGhostSliceOn(side, {0});
	loop_over_interior_indexes<3>(local_ghosts, [&](const std::array<int, 3> &coord) { local_ghosts[coord] -= 1.0 / 3.0 * local_slice[coord]; });
}
/**
 * @brief Fill ghost cells for a finer neighbor
 *
 * This is only part of the interpolation, values are also needed from the neighbor patch.
 * These values are added in FillGhostCaellsForLocalWithCoarseNbr
 *
 * @param local_view the local patch data
 * @param nbr_view the neighbor patch adata
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of this patches side that the neighbor lies on
 */
void FillGhostCellsForFineNbr(const PatchView<const double, 3> &local_view,
                              const PatchView<const double, 3> &nbr_view,
                              Side<3>                           side,
                              Orthant<2>                        orthant)
{
	auto [offset_i, offset_j]         = getOffset(local_view.getEnd(), side, orthant);
	View<const double, 3> local_slice = local_view.getSliceOn(side, {0});
	View<double, 3>       nbr_ghosts  = nbr_view.getGhostSliceOn(side.opposite(), {0});

	for (int c = local_slice.getStart()[2]; c <= local_slice.getEnd()[2]; c++) {
		for (int j = local_slice.getStart()[1]; j <= local_slice.getEnd()[1]; j++) {
			for (int i = local_slice.getStart()[0]; i <= local_slice.getEnd()[0]; i++) {
				nbr_ghosts(i, j, c) += 4.0 / 6.0 * local_slice((i + offset_i) / 2, (j + offset_j) / 2, c);
			}
		}
	}
}
/**
 * @brief Add in the values need for this patch's ghost cells when there is a coarse neighbor
 *
 * These values compliment FIllGhostCellsForFineNbr
 *
 * @param pinfo the patchinfo object
 * @param view the local patch view
 * @param side the side of the patch that the neighbor is on
 */
void FillGhostCellsForLocalWithCoarseNbr(const PatchInfo<3> &pinfo, PatchView<const double, 3> view, Side<3> side)
{
	View<const double, 3>   local_slice   = view.getSliceOn(side, {0});
	View<double, 3>         local_ghosts  = view.getGhostSliceOn(side, {0});
	const CoarseNbrInfo<2> &nbr_info      = pinfo.getCoarseNbrInfo(side);
	auto [coarse_start_i, coarse_start_j] = getOffset(view.getEnd(), side, nbr_info.orth_on_coarse);

	for (int c = local_slice.getStart()[2]; c <= local_slice.getEnd()[2]; c++) {
		for (int j = local_slice.getStart()[1]; j <= local_slice.getEnd()[1]; j++) {
			int offset_j = (j + coarse_start_j) % 2 == 0 ? (j + 1) : (j - 1);
			for (int i = local_slice.getStart()[0]; i <= local_slice.getEnd()[0]; i++) {
				local_ghosts(i, j, c) += 5.0 / 6.0 * local_slice(i, j, c);
				int offset_i = (i + coarse_start_i) % 2 == 0 ? (i + 1) : (i - 1);
				local_ghosts(i, j, c) -= 1.0 / 6.0 * local_slice(offset_i, j, c);
				local_ghosts(i, j, c) -= 1.0 / 6.0 * local_slice(i, offset_j, c);
				local_ghosts(i, j, c) -= 1.0 / 6.0 * local_slice(offset_i, offset_j, c);
			}
		}
	}
}

////////
// EDGES
////////

/**
 * @brief Simple copy of values
 *
 * @param local_view the neighbor data
 * @param nbr_view  the local data
 * @param edge the edge that the neighbor is on
 */
void FillGhostCellsForNormalEdgeNbr(const PatchView<const double, 3> &local_view, const PatchView<const double, 3> &nbr_view, Edge edge)
{
	View<const double, 2> local_slice = local_view.getSliceOn(edge, {0, 0});
	View<double, 2>       nbr_ghosts  = nbr_view.getGhostSliceOn(edge.opposite(), {0, 0});
	loop_over_interior_indexes<2>(nbr_ghosts, [&](const std::array<int, 2> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
}
/**
 * @brief Fill ghost cells for coarser neighbor
 *
 * This is only part of the interpolation, values are also needed from the neighbor patch.
 * These values are added in FillGhostCaellsForLocalWithFineEdgeNbr
 *
 * @param local_view local patch data
 * @param nbr_view neighbor patch data
 * @param edge the edge of the patch that the neighbor patch is on
 * @param orthant the orthant of the neighbors edge that this patch lies on
 */
void FillGhostCellsForCoarseEdgeNbr(const PatchView<const double, 3> &local_view,
                                    const PatchView<const double, 3> &nbr_view,
                                    Edge                              edge,
                                    Orthant<1>                        orthant)
{
	View<const double, 2> local_slice = local_view.getSliceOn(edge, {0, 0});
	View<double, 2>       nbr_ghosts  = nbr_view.getGhostSliceOn(edge.opposite(), {0, 0});
	int                   offset      = 0;
	if (orthant == Orthant<1>::upper()) {
		offset = local_slice.getEnd()[0] + 1;
	}
	for (int c = local_slice.getStart()[1]; c <= local_slice.getEnd()[1]; c++) {
		for (int i = local_slice.getStart()[0]; i <= local_slice.getEnd()[0]; i++) {
			nbr_ghosts((i + offset) / 2, c) += 2.0 / 3.0 * local_slice(i, c);
		}
	}
}
/**
 * @brief Fill ghost cells for a finer neighbor
 *
 * This is only part of the interpolation, values are also needed from the neighbor patch.
 * These values are added in FillGhostCaellsForLocalWithCoarseEdgeNbr
 *
 * @param local_view the local patch data
 * @param nbr_view the neighbor patch adata
 * @param edge the edge that the neighbor patch is on
 * @param orthant the orthant of this patches edge that the neighbor lies on
 */
void FillGhostCellsForFineEdgeNbr(const PatchView<const double, 3> &local_view,
                                  const PatchView<const double, 3> &nbr_view,
                                  Edge                              edge,
                                  Orthant<1>                        orthant)
{
	View<const double, 2> local_slice = local_view.getSliceOn(edge, {0, 0});
	View<double, 2>       nbr_ghosts  = nbr_view.getGhostSliceOn(edge.opposite(), {0, 0});
	int                   offset      = 0;
	if (orthant == Orthant<1>::upper()) {
		offset = local_slice.getEnd()[0] + 1;
	}
	for (int c = local_slice.getStart()[1]; c <= local_slice.getEnd()[1]; c++) {
		for (int i = local_slice.getStart()[0]; i <= local_slice.getEnd()[0]; i++) {
			nbr_ghosts(i, c) += 2.0 / 3.0 * local_slice((i + offset) / 2, c);
		}
	}
}
/**
 * @brief Add in the values need for this patch's ghost cells when there is a coarse neighbor
 *
 * These values compliment FIllGhostCellsForFineEdgeNbr
 *
 * @param pinfo the patchinfo object
 * @param view the local patch view
 * @param edge the edge of the patch that the neighbor is on
 */
void FillGhostCellsForLocalWithCoarseEdgeNbr(const PatchInfo<3> &pinfo, const PatchView<const double, 3> &view, Edge edge)
{
	View<const double, 2> local_slice  = view.getSliceOn(edge, {0, 0});
	View<double, 2>       local_ghosts = view.getGhostSliceOn(edge, {0, 0});

	int offset = 0;
	if (pinfo.getCoarseNbrInfo(edge).orth_on_coarse == Orthant<1>::upper()) {
		offset = local_slice.getEnd()[0] + 1;
	}
	for (int c = local_slice.getStart()[1]; c <= local_slice.getEnd()[1]; c++) {
		for (int i = local_slice.getStart()[0]; i <= local_slice.getEnd()[0]; i++) {
			local_ghosts(i, c) += 2.0 / 3.0 * local_slice(i, c);

			if ((i + offset) % 2 == 0) {
				local_ghosts(i + 1, c) += -1.0 / 3.0 * local_slice(i, c);
			} else {
				local_ghosts(i - 1, c) += -1.0 / 3.0 * local_slice(i, c);
			}
		}
	}
}
/**
 * @brief Add in the values need for this patches ghost cells when there is a fine neighbor
 *
 * These values compliment FIllGhostCellsForCoarseEdgeNbr
 *
 * @param local_data the patch data
 * @param edge the edge of the patch that the neighbor is on
 */
void FillGhostCellsForLocalWithFineEdgeNbr(const PatchView<const double, 3> &view, Edge edge)
{
	View<const double, 2> local_slice  = view.getSliceOn(edge, {0, 0});
	View<double, 2>       local_ghosts = view.getGhostSliceOn(edge, {0, 0});
	loop_over_interior_indexes<2>(local_ghosts, [&](const std::array<int, 2> &coord) { local_ghosts[coord] += -1.0 / 3.0 * local_slice[coord]; });
}

//////////
// CORNERS
//////////

/**
 * @brief Simple copy of values
 *
 * @param local_view the neighbor data
 * @param nbr_view  the local data
 * @param corner the corner that the neighbor is on
 */
void FillGhostCellsForNormalCornerNbr(const PatchView<const double, 3> &local_view, const PatchView<const double, 3> &nbr_view, Corner<3> corner)
{
	View<const double, 1> local_slice = local_view.getSliceOn(corner, {0, 0, 0});
	View<double, 1>       nbr_ghost   = nbr_view.getGhostSliceOn(corner.opposite(), {0, 0, 0});
	loop_over_interior_indexes<1>(local_slice, [&](const std::array<int, 1> &coord) { nbr_ghost[coord] = local_slice[coord]; });
}

/**
 * @brief Fill ghost cells for coarser neighbor
 *
 * This is only part of the interpolation, values are also needed from the neighbor patch.
 * These values are added in FillGhostCaellsForLocalWithFineCornerNbr
 *
 * @param local_view local patch data
 * @param nbr_view neighbor patch data
 * @param corner the corner of the patch that the neighbor patch is on
 */
void FillGhostCellsForCoarseCornerNbr(const PatchView<const double, 3> &local_view, const PatchView<const double, 3> &nbr_view, Corner<3> corner)
{
	View<const double, 1> local_slice = local_view.getSliceOn(corner, {0, 0, 0});
	View<double, 1>       nbr_ghost   = nbr_view.getGhostSliceOn(corner.opposite(), {0, 0, 0});
	loop_over_interior_indexes<1>(local_slice, [&](const std::array<int, 1> &coord) { nbr_ghost[coord] += 4.0 * local_slice[coord] / 3.0; });
}
/**
 * @brief Fill ghost cells for a finer neighbor
 *
 * This is only part of the interpolation, values are also needed from the neighbor patch.
 * These values are added in FillGhostCaellsForLocalWithCoarseEdgeNbr
 *
 * @param local_view the local patch data
 * @param nbr_view the neighbor patch adata
 * @param edge the edge that the neighbor patch is on
 * @param orthant the orthant of this patches edge that the neighbor lies on
 */
void FillGhostCellsForFineCornerNbr(const PatchView<const double, 3> &local_view, const PatchView<const double, 3> &nbr_view, Corner<3> corner)
{
	View<const double, 1> local_slice = local_view.getSliceOn(corner, {0, 0, 0});
	View<double, 1>       nbr_ghost   = nbr_view.getGhostSliceOn(corner.opposite(), {0, 0, 0});
	loop_over_interior_indexes<1>(local_slice, [&](const std::array<int, 1> &coord) { nbr_ghost[coord] += 2.0 * local_slice[coord] / 3.0; });
}
/**
 * @brief Add in the values need for this patch's ghost cells when there is a coarse neighbor
 *
 * These values compliment FIllGhostCellsForFineCornerNbr
 *
 * @param view the local patch view
 * @param corner the corner of the patch that the neighbor is on
 */
void FillGhostCellsForLocalWithCoarseCornerNbr(const PatchView<const double, 3> &view, Corner<3> corner)
{
	View<const double, 1> local_slice = view.getSliceOn(corner, {0, 0, 0});
	View<double, 1>       local_ghost = view.getGhostSliceOn(corner, {0, 0, 0});
	loop_over_interior_indexes<1>(local_slice, [&](const std::array<int, 1> &coord) { local_ghost[coord] += local_slice[coord] / 3.0; });
}
/**
 * @brief Add in the values need for this patches ghost cells when there is a fine neighbor
 *
 * These values compliment FIllGhostCellsForCoarseCornerNbr
 *
 * @param view the patch view
 * @param corner the corner of the patch that the neighbor is on
 */
void FillGhostCellsForLocalWithFineCornerNbr(const PatchView<const double, 3> &view, Corner<3> corner)
{
	View<const double, 1> local_slice = view.getSliceOn(corner, {0, 0, 0});
	View<double, 1>       local_ghost = view.getGhostSliceOn(corner, {0, 0, 0});
	loop_over_interior_indexes<1>(local_slice, [&](const std::array<int, 1> &coord) { local_ghost[coord] += -local_slice[coord] / 3.0; });
}
/**
 * @brief Add in extra information needed from the local patch on the sides
 *
 * @param pinfo the pinfo object
 * @param view the patch view
 */
void FillLocalGhostCellsOnSides(const PatchInfo<3> &pinfo, const PatchView<const double, 3> &view)
{
	for (Side<3> side : Side<3>::getValues()) {
		if (pinfo.hasNbr(side)) {
			switch (pinfo.getNbrType(side)) {
				case NbrType::Normal:
					// do nothing
					break;
				case NbrType::Coarse:
					FillGhostCellsForLocalWithCoarseNbr(pinfo, view, side);
					break;
				case NbrType::Fine:
					FillGhostCellsForLocalWithFineNbr(view, side);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
/**
 * @brief Add in extra information needed from the local patch on the edges
 *
 * @param pinfo the pinfo object
 * @param view the patch view
 */
void FillLocalGhostCellsOnEdges(const PatchInfo<3> &pinfo, const PatchView<const double, 3> &view)
{
	for (Edge edge : Edge::getValues()) {
		if (pinfo.hasNbr(edge)) {
			switch (pinfo.getNbrType(edge)) {
				case NbrType::Normal:
					// do nothing
					break;
				case NbrType::Coarse:
					FillGhostCellsForLocalWithCoarseEdgeNbr(pinfo, view, edge);
					break;
				case NbrType::Fine:
					FillGhostCellsForLocalWithFineEdgeNbr(view, edge);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
/**
 * @brief Add in extra information needed from the local patch on the corners
 *
 * @param pinfo the pinfo object
 * @param view the patch view
 */
void FillLocalGhostCellsOnCorners(const PatchInfo<3> &pinfo, const PatchView<const double, 3> &view)
{
	for (Corner<3> corner : Corner<3>::getValues()) {
		if (pinfo.hasNbr(corner)) {
			switch (pinfo.getNbrType(corner)) {
				case NbrType::Normal:
					// do nothing
					break;
				case NbrType::Coarse:
					FillGhostCellsForLocalWithCoarseCornerNbr(view, corner);
					break;
				case NbrType::Fine:
					FillGhostCellsForLocalWithFineCornerNbr(view, corner);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
} // namespace

void TriLinearGhostFiller::fillGhostCellsForNbrPatch(const PatchInfo<3> &              pinfo,
                                                     const PatchView<const double, 3> &local_view,
                                                     const PatchView<const double, 3> &nbr_view,
                                                     Side<3>                           side,
                                                     NbrType                           nbr_type,
                                                     Orthant<2>                        orthant) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostCellsForNormalNbr(local_view, nbr_view, side);
			break;
		case NbrType::Coarse:
			FillGhostCellsForCoarseNbr(local_view, nbr_view, side, orthant);
			break;
		case NbrType::Fine:
			FillGhostCellsForFineNbr(local_view, nbr_view, side, orthant);
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void TriLinearGhostFiller::fillGhostCellsForEdgeNbrPatch(const PatchInfo<3> &              pinfo,
                                                         const PatchView<const double, 3> &local_view,
                                                         const PatchView<const double, 3> &nbr_view,
                                                         Edge                              edge,
                                                         NbrType                           nbr_type,
                                                         Orthant<1>                        orthant_on_coarse) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostCellsForNormalEdgeNbr(local_view, nbr_view, edge);
			break;
		case NbrType::Coarse:
			FillGhostCellsForCoarseEdgeNbr(local_view, nbr_view, edge, orthant_on_coarse);
			break;
		case NbrType::Fine:
			FillGhostCellsForFineEdgeNbr(local_view, nbr_view, edge, orthant_on_coarse);
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void TriLinearGhostFiller::fillGhostCellsForCornerNbrPatch(const PatchInfo<3> &              pinfo,
                                                           const PatchView<const double, 3> &local_view,
                                                           const PatchView<const double, 3> &nbr_view,
                                                           Corner<3>                         corner,
                                                           NbrType                           nbr_type) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostCellsForNormalCornerNbr(local_view, nbr_view, corner);
			break;
		case NbrType::Coarse:
			FillGhostCellsForCoarseCornerNbr(local_view, nbr_view, corner);
			break;
		case NbrType::Fine:
			FillGhostCellsForFineCornerNbr(local_view, nbr_view, corner);
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void TriLinearGhostFiller::fillGhostCellsForLocalPatch(const PatchInfo<3> &pinfo, const PatchView<const double, 3> &local_view) const
{
	switch (this->getFillType()) {
		case GhostFillingType::Corners: // Fill corners, edges, and faces
			FillLocalGhostCellsOnCorners(pinfo, local_view);
			[[fallthrough]];
		case GhostFillingType::Edges: // Fill edges and faces
			FillLocalGhostCellsOnEdges(pinfo, local_view);
			[[fallthrough]];
		case GhostFillingType::Faces: // Fill faces
			FillLocalGhostCellsOnSides(pinfo, local_view);
			break;
		default:
			throw RuntimeError("Unsupported GhostFillingType");
	}
}
TriLinearGhostFiller::TriLinearGhostFiller(std::shared_ptr<const Domain<3>> domain, GhostFillingType fill_type) : MPIGhostFiller<3>(domain, fill_type)
{
	for (int n : domain->getNs()) {
		if (n % 2 != 0) {
			throw RuntimeError("TriLinearGhostFiller only supports even number patch sizes");
		}
	}
}

} // namespace ThunderEgg