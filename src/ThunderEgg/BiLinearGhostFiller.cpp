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

#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/RuntimeError.h>
namespace ThunderEgg
{
BiLinearGhostFiller::BiLinearGhostFiller(std::shared_ptr<const Domain<2>> domain, GhostFillingType fill_type) : MPIGhostFiller<2>(domain, fill_type)
{
}
namespace
{
/**
 * @brief This is just a simple copy of values
 *
 * @param local_datas the local pach
 * @param nbr_datas the neighbor patch
 * @param side the side that the neighbor patch is on
 */
void FillGhostForNormalNbr(const std::vector<ComponentView<2>> &local_datas, const std::vector<ComponentView<2>> &nbr_datas, const Side<2> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		ConstView<1> local_slice = local_datas[c].getSliceOn(side, {0});
		View<1>      nbr_ghosts  = nbr_datas[c].getGhostSliceOn(side.opposite(), {0});
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 1> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
	}
}
/**
 * @brief Fill the ghost values for a coarse neighbor
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be needed.
 *
 * Those values are added in FillLocalGhostsForFineNbr
 *
 * @param local_datas the local patch
 * @param nbr_datas the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the coarser neighbor face that this patch lies on
 */
void FillGhostForCoarseNbr(const std::vector<ComponentView<2>> &local_datas,
                           const std::vector<ComponentView<2>> &nbr_datas,
                           const Side<2>                        side,
                           const Orthant<1>                     orthant)
{
	int offset = 0;
	if (orthant == Orthant<1>::upper()) {
		offset = local_datas[0].getLengths()[!side.getAxisIndex()];
	}
	for (size_t c = 0; c < local_datas.size(); c++) {
		ConstView<1> local_slice = local_datas[c].getSliceOn(side, {0});
		View<1>      nbr_ghosts  = nbr_datas[c].getGhostSliceOn(side.opposite(), {0});
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 1> &coord) {
			nbr_ghosts[{(coord[0] + offset) / 2}] += 2.0 / 3.0 * local_slice[coord];
		});
	}
}
/**
 * @brief Fill the ghost values for a fine neighbor
 *
 * This is only part of the value. Values from the interior of the fine neighbor will also be needed.
 *
 * Those values are added in FillLocalGhostsForCoarseNbr
 *
 * @param local_datas the local patch
 * @param nbr_datas the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the this patches face that the finer neighbor patch lies on
 */
void FillGhostForFineNbr(const std::vector<ComponentView<2>> &local_datas,
                         const std::vector<ComponentView<2>> &nbr_datas,
                         const Side<2>                        side,
                         const Orthant<1>                     orthant)
{
	int offset = 0;
	if (orthant == Orthant<1>::upper()) {
		offset = local_datas[0].getLengths()[!side.getAxisIndex()];
	}
	for (size_t c = 0; c < local_datas.size(); c++) {
		ConstView<1> local_slice = local_datas[c].getSliceOn(side, {0});
		View<1>      nbr_ghosts  = nbr_datas[c].getGhostSliceOn(side.opposite(), {0});
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 1> &coord) {
			nbr_ghosts[coord] += 2.0 / 3.0 * local_slice[{(coord[0] + offset) / 2}];
		});
	}
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a coarser neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForFineNbr
 *
 * @param pinfo the patchinfo object
 * @param local_data the patch data
 * @param side the side that the neighbor patch is on
 */
void FillLocalGhostsForCoarseNbr(const PatchInfo<2> &pinfo, const ComponentView<2> &local_data, const Side<2> side)
{
	ConstView<1> local_slice  = local_data.getSliceOn(side, {0});
	View<1>      local_ghosts = local_data.getGhostSliceOn(side, {0});
	int          offset       = 0;
	if (pinfo.getCoarseNbrInfo(side).orth_on_coarse == Orthant<1>::upper()) {
		offset = local_data.getLengths()[!side.getAxisIndex()];
	}
	nested_loop<1>(local_ghosts.getStart(), local_ghosts.getEnd(), [&](const std::array<int, 1> &coord) {
		local_ghosts[coord] += 2.0 / 3.0 * local_slice[coord];
		if ((coord[0] + offset) % 2 == 0) {
			local_ghosts[{coord[0] + 1}] += -1.0 / 3.0 * local_slice[coord];
		} else {
			local_ghosts[{coord[0] - 1}] += -1.0 / 3.0 * local_slice[coord];
		}
	});
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a finer neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForCoarseNbr
 *
 * @param local_data the patch data
 * @param side the side that the neighbor patch is on
 */
void FillLocalGhostsForFineNbr(const ComponentView<2> &local_data, const Side<2> side)
{
	ConstView<1> local_slice  = local_data.getSliceOn(side, {0});
	View<1>      local_ghosts = local_data.getGhostSliceOn(side, {0});
	nested_loop<1>(
	local_ghosts.getStart(), local_ghosts.getEnd(), [&](const std::array<int, 1> &coord) { local_ghosts[coord] += -1.0 / 3.0 * local_slice[coord]; });
}
/**
 * @brief This is just a simple copy of values
 *
 * @param local_datas the local pach
 * @param nbr_datas the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForCornerNormalNbr(const std::vector<ComponentView<2>> &local_datas, const std::vector<ComponentView<2>> &nbr_datas, Corner<2> corner)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		ConstView<0> local_slice = local_datas[c].getSliceOn(corner, {0, 0});
		View<0>      nbr_ghost   = nbr_datas[c].getGhostSliceOn(corner.opposite(), {0, 0});
		nbr_ghost[{}]            = local_slice[{}];
	}
}
/**
 * @brief Fill the ghost values for a coarse neighbor
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be needed.
 *
 * Those values are added in FillLocalGhostsForCornerFineNbr
 *
 * @param local_datas the local patch
 * @param nbr_datas the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForCornerCoarseNbr(const std::vector<ComponentView<2>> &local_datas, const std::vector<ComponentView<2>> &nbr_datas, Corner<2> corner)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		View<0> nbr_ghosts = nbr_datas[c].getGhostSliceOn(corner.opposite(), {0, 0});
		nbr_ghosts[{}] += 4.0 * local_datas[c].getSliceOn(corner, {0, 0})[{}] / 3.0;
	}
}
/**
 * @brief Fill the ghost values for a fine neighbor
 *
 * This is only part of the value. Values from the interior of the fine neighbor will also be needed.
 *
 * Those values are added in FillLocalGhostsForCornerCoarseNbr
 *
 * @param local_datas the local patch
 * @param nbr_datas the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void FillGhostForCornerFineNbr(const std::vector<ComponentView<2>> &local_datas, const std::vector<ComponentView<2>> &nbr_datas, Corner<2> corner)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		ConstView<0> local_slice = local_datas[c].getSliceOn(corner, {0, 0});
		View<0>      nbr_ghosts  = nbr_datas[c].getGhostSliceOn(corner.opposite(), {0, 0});
		nbr_ghosts[{}] += 2.0 * local_slice[{}] / 3.0;
	}
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a coarser neighbor
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be needed.
 *
 * Those values compliment FillGhostForCornerFineNbr
 *
 * @param local_data the patch data
 * @param corner the corner that the neighbor patch is on
 */
void FillLocalGhostsForCornerCoarseNbr(const ComponentView<2> &local_data, Corner<2> corner)
{
	ConstView<0> local_slice  = local_data.getSliceOn(corner, {0, 0});
	View<0>      local_ghosts = local_data.getGhostSliceOn(corner, {0, 0});
	local_ghosts[{}] += local_slice[{}] / 3.0;
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a finer neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForCornerCoarseNbr
 *
 * @param local_data the patch data
 * @param corner the corner that the neighbor patch is on
 */
void FillLocalGhostsForCornerFineNbr(const ComponentView<2> &local_data, Corner<2> corner)
{
	ConstView<0> local_slice  = local_data.getSliceOn(corner, {0, 0});
	View<0>      local_ghosts = local_data.getGhostSliceOn(corner, {0, 0});
	local_ghosts[{}] += -local_slice[{}] / 3.0;
}
/**
 * @brief Add in extra information needed from the local patch on the sides
 *
 * @param pinfo the pinfo object
 * @param local_data the patch data
 */
void FillLocalGhostCellsOnSides(const PatchInfo<2> &pinfo, const ComponentView<2> &local_data)
{
	for (Side<2> side : Side<2>::getValues()) {
		if (pinfo.hasNbr(side)) {
			switch (pinfo.getNbrType(side)) {
				case NbrType::Normal:
					// nothing needs to be done
					break;
				case NbrType::Coarse:
					FillLocalGhostsForCoarseNbr(pinfo, local_data, side);
					break;
				case NbrType::Fine:
					FillLocalGhostsForFineNbr(local_data, side);
					break;
				default:
					throw RuntimeError("Unsupported Nbr Type");
			}
		}
	}
}
/**
 * @brief Add in extra information needed from the local patch on the corner
 *
 * @param pinfo the pinfo object
 * @param local_data the patch data
 */
void FillLocalGhostCellsOnCorners(const PatchInfo<2> &pinfo, const ComponentView<2> &local_data)
{
	for (Corner<2> corner : Corner<2>::getValues()) {
		if (pinfo.hasNbr(corner)) {
			switch (pinfo.getNbrType(corner)) {
				case NbrType::Normal:
					// nothing needs to be done
					break;
				case NbrType::Coarse:
					FillLocalGhostsForCornerCoarseNbr(local_data, corner);
					break;
				case NbrType::Fine:
					FillLocalGhostsForCornerFineNbr(local_data, corner);
					break;
				default:
					throw RuntimeError("Unsupported Nbr Type");
			}
		}
	}
}
} // namespace
void BiLinearGhostFiller::fillGhostCellsForNbrPatch(const PatchInfo<2> &                 pinfo,
                                                    const std::vector<ComponentView<2>> &local_datas,
                                                    std::vector<ComponentView<2>> &      nbr_datas,
                                                    Side<2>                              side,
                                                    NbrType                              nbr_type,
                                                    Orthant<1>                           orthant_on_coarse) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostForNormalNbr(local_datas, nbr_datas, side);
			break;
		case NbrType::Coarse:
			FillGhostForCoarseNbr(local_datas, nbr_datas, side, orthant_on_coarse);
			break;
		case NbrType::Fine:
			FillGhostForFineNbr(local_datas, nbr_datas, side, orthant_on_coarse);
			break;
		default:
			throw RuntimeError("Unsupported Nbr Type");
	}
}
void BiLinearGhostFiller::fillGhostCellsForEdgeNbrPatch(const PatchInfo<2> &                 pinfo,
                                                        const std::vector<ComponentView<2>> &local_datas,
                                                        std::vector<ComponentView<2>> &      nbr_datas,
                                                        Edge                                 edge,
                                                        NbrType                              nbr_type,
                                                        Orthant<1>                           orthant_on_coarse) const
{
	// 2D, edges not needed
}

void BiLinearGhostFiller::fillGhostCellsForCornerNbrPatch(const PatchInfo<2> &                 pinfo,
                                                          const std::vector<ComponentView<2>> &local_datas,
                                                          std::vector<ComponentView<2>> &      nbr_datas,
                                                          Corner<2>                            corner,
                                                          NbrType                              nbr_type) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostForCornerNormalNbr(local_datas, nbr_datas, corner);
			break;
		case NbrType::Coarse:
			FillGhostForCornerCoarseNbr(local_datas, nbr_datas, corner);
			break;
		case NbrType::Fine:
			FillGhostForCornerFineNbr(local_datas, nbr_datas, corner);
			break;
		default:
			throw RuntimeError("Unsupported Nbr Type");
	}
}
void BiLinearGhostFiller::fillGhostCellsForLocalPatch(const PatchInfo<2> &pinfo, std::vector<ComponentView<2>> &local_datas) const
{
	for (const ComponentView<2> &local_data : local_datas) {
		switch (this->getFillType()) {
			case GhostFillingType::Corners: // Fill corners and faces
				FillLocalGhostCellsOnCorners(pinfo, local_data);
				[[fallthrough]];
			case GhostFillingType::Faces:
				FillLocalGhostCellsOnSides(pinfo, local_data);
				break;
			default:
				throw RuntimeError("Unsupported GhostFillingType");
		}
	}
}
} // namespace ThunderEgg