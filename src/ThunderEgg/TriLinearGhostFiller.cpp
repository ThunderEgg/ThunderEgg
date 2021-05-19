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
 * @param ns the dimensions of the patch
 * @param orth the orthant
 * @return std::array<int,2> the offsets
 */
std::array<int, 2> getOffset(const std::array<int, 3> ns, Side<3> s, Orthant<2> orth)
{
	std::array<int, 2> offset = {0, 0};
	for (size_t i = 0; i < s.getAxisIndex(); i++) {
		if (orth.isHigherOnAxis(i)) {
			offset[i] = ns[i];
		}
	}
	for (size_t i = s.getAxisIndex() + 1; i < 3; i++) {
		if (orth.isHigherOnAxis(i - 1)) {
			offset[i - 1] = ns[i];
		}
	}
	return offset;
}
void FillGhostCellsForNormalNbr(const std::vector<LocalData<3>> &local_datas, std::vector<LocalData<3>> &nbr_datas, Side<3> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		LocalData<2> local_slice = local_datas[c].getSliceOn(side, {0});
		LocalData<2> nbr_ghosts  = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		nested_loop<2>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 2> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
	}
}
void FillGhostCellsForCoarseNbr(const PatchInfo<3> &             pinfo,
                                const std::vector<LocalData<3>> &local_datas,
                                std::vector<LocalData<3>> &      nbr_datas,
                                Side<3>                          side,
                                Orthant<2>                       orthant)
{
	std::array<int, 2> offset = getOffset(pinfo.ns, side, orthant);
	for (size_t c = 0; c < local_datas.size(); c++) {
		LocalData<2> local_slice = local_datas[c].getSliceOn(side, {0});
		LocalData<2> nbr_ghosts  = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		nested_loop<2>(local_slice.getStart(), local_slice.getEnd(), [&](const std::array<int, 2> &coord) {
			std::array<int, 2> coarse_coord;
			for (int i = 0; i < 2; i++) {
				coarse_coord[i] = (coord[i] + offset[i]) / 2;
			}
			nbr_ghosts[coarse_coord] += 1.0 / 3.0 * local_slice[coord];
		});
	}
}
void FillGhostCellsForLocalWithFineNbr(LocalData<3> local_data, Side<3> side)
{
	auto local_slice  = local_data.getSliceOn(side, {0});
	auto local_ghosts = local_data.getSliceOn(side, {-1});
	nested_loop<2>(
	local_ghosts.getStart(), local_ghosts.getEnd(), [&](const std::array<int, 2> &coord) { local_ghosts[coord] -= 1.0 / 3.0 * local_slice[coord]; });
}
void FillGhostCellsForFineNbr(const PatchInfo<3> &             pinfo,
                              const std::vector<LocalData<3>> &local_datas,
                              std::vector<LocalData<3>> &      nbr_datas,
                              Side<3>                          side,
                              Orthant<2>                       orthant)
{
	auto               nbr_info = pinfo.getFineNbrInfo(side);
	std::array<int, 2> offset   = getOffset(pinfo.ns, side, orthant);
	for (size_t c = 0; c < local_datas.size(); c++) {
		LocalData<2> local_slice = local_datas[c].getSliceOn(side, {0});
		LocalData<2> nbr_ghosts  = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		nested_loop<2>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 2> &coord) {
			std::array<int, 2> coarse_coord;
			for (int i = 0; i < 2; i++) {
				coarse_coord[i] = (coord[i] + offset[i]) / 2;
			}
			nbr_ghosts[coord] += 4.0 / 6.0 * local_slice[coarse_coord];
		});
	}
}
void FillGhostCellsForLocalWithCoarseNbr(const PatchInfo<3> &pinfo, LocalData<3> local_data, Side<3> side)
{
	LocalData<2>       local_slice  = local_data.getSliceOn(side, {0});
	LocalData<2>       local_ghosts = local_data.getSliceOn(side, {-1});
	auto               nbr_info     = pinfo.getCoarseNbrInfo(side);
	std::array<int, 2> offset       = getOffset(pinfo.ns, side, nbr_info.orth_on_coarse);
	nested_loop<2>(local_ghosts.getStart(), local_ghosts.getEnd(), [&](const std::array<int, 2> &coord) {
		std::array<int, 2> offset_coord;
		for (int i = 0; i < 2; i++) {
			if ((coord[i] + offset[i]) % 2 == 0) {
				offset_coord[i] = coord[i] + 1;
			} else {
				offset_coord[i] = coord[i] - 1;
			}
		}
		local_ghosts[coord] += 5.0 / 6.0 * local_slice[coord];
		local_ghosts[coord] -= 1.0 / 6.0 * local_slice[{offset_coord[0], coord[1]}];
		local_ghosts[coord] -= 1.0 / 6.0 * local_slice[{coord[0], offset_coord[1]}];
		local_ghosts[coord] -= 1.0 / 6.0 * local_slice[offset_coord];
	});
}
void FillGhostCellsForNormalEdgeNbr(const std::vector<LocalData<3>> &local_datas, const std::vector<LocalData<3>> &nbr_datas, Edge edge)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		LocalData<1> local_slice = local_datas[c].getSliceOn(edge, {0, 0});
		LocalData<1> nbr_ghosts  = nbr_datas[c].getSliceOn(edge.opposite(), {-1, -1});
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 1> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
	}
}
void FillGhostCellsForCoarseEdgeNbr(const std::vector<LocalData<3>> &local_datas,
                                    const std::vector<LocalData<3>> &nbr_datas,
                                    Edge                             edge,
                                    Orthant<1>                       orthant)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		LocalData<1> local_slice = local_datas[c].getSliceOn(edge, {0, 0});
		LocalData<1> nbr_ghosts  = nbr_datas[c].getSliceOn(edge.opposite(), {-1, -1});
		int          offset      = 0;
		if (orthant == Orthant<1>::upper()) {
			offset = local_slice.getLengths()[0];
		}
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 1> &coord) {
			nbr_ghosts[{(coord[0] + offset) / 2}] += 2.0 / 3.0 * local_slice[coord];
		});
	}
}
void FillGhostCellsForFineEdgeNbr(const std::vector<LocalData<3>> &local_datas,
                                  const std::vector<LocalData<3>> &nbr_datas,
                                  Edge                             edge,
                                  Orthant<1>                       orthant)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		LocalData<1> local_slice = local_datas[c].getSliceOn(edge, {0, 0});
		LocalData<1> nbr_ghosts  = nbr_datas[c].getSliceOn(edge.opposite(), {-1, -1});
		int          offset      = 0;
		if (orthant == Orthant<1>::upper()) {
			offset = local_slice.getLengths()[0];
		}
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 1> &coord) {
			nbr_ghosts[coord] += 2.0 / 3.0 * local_slice[{(coord[0] + offset) / 2}];
		});
	}
}
void FillGhostCellsForLocalWithCoarseEdgeNbr(const PatchInfo<3> &pinfo, const LocalData<3> &local_data, Edge side)
{
	LocalData<1> local_slice  = local_data.getSliceOn(side, {0, 0});
	LocalData<1> local_ghosts = local_data.getSliceOn(side, {-1, -1});
	nested_loop<1>(local_ghosts.getStart(), local_ghosts.getEnd(), [&](const std::array<int, 1> &coord) {
		local_ghosts[coord] += 2.0 / 3.0 * local_slice[coord];

		int offset = 0;
		if (pinfo.getCoarseNbrInfo(side).orth_on_coarse == Orthant<1>::upper()) {
			offset = local_slice.getLengths()[0];
		}

		if ((coord[0] + offset) % 2 == 0) {
			local_ghosts[{coord[0] + 1}] += -1.0 / 3.0 * local_slice[coord];
		} else {
			local_ghosts[{coord[0] - 1}] += -1.0 / 3.0 * local_slice[coord];
		}
	});
}
void FillGhostCellsForLocalWithFineEdgeNbr(const LocalData<3> &local_data, Edge side)
{
	LocalData<1> local_slice  = local_data.getSliceOn(side, {0, 0});
	LocalData<1> local_ghosts = local_data.getSliceOn(side, {-1, -1});
	nested_loop<1>(
	local_ghosts.getStart(), local_ghosts.getEnd(), [&](const std::array<int, 1> &coord) { local_ghosts[coord] += -1.0 / 3.0 * local_slice[coord]; });
}
void FillGhostCellsForNormalCornerNbr(const std::vector<LocalData<3>> &local_datas, const std::vector<LocalData<3>> &nbr_datas, Corner<3> corner)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto local_slice = local_datas[c].getSliceOn(corner, {0, 0, 0});
		auto nbr_ghost   = nbr_datas[c].getSliceOn(corner.opposite(), {-1, -1, -1});
		nbr_ghost[{}]    = local_slice[{}];
	}
}
void FillGhostCellsForCoarseCornerNbr(const std::vector<LocalData<3>> &local_datas, const std::vector<LocalData<3>> &nbr_datas, Corner<3> corner)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto nbr_ghosts = nbr_datas[c].getSliceOn(corner.opposite(), {-1, -1, -1});
		nbr_ghosts[{}] += 4.0 * local_datas[c].getSliceOn(corner, {0, 0, 0})[{}] / 3.0;
	}
}
void FillGhostCellsForFineCornerNbr(const std::vector<LocalData<3>> &local_datas, const std::vector<LocalData<3>> &nbr_datas, Corner<3> corner)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto local_slice = local_datas[c].getSliceOn(corner, {0, 0, 0});
		auto nbr_ghosts  = nbr_datas[c].getSliceOn(corner.opposite(), {-1, -1, -1});
		nbr_ghosts[{}] += 2.0 * local_slice[{}] / 3.0;
	}
}
void FillGhostCellsForLocalWithCoarseCornerNbr(const PatchInfo<3> &pinfo, const LocalData<3> &local_data, Corner<3> corner)
{
	auto local_slice  = local_data.getSliceOn(corner, {0, 0, 0});
	auto local_ghosts = local_data.getSliceOn(corner, {-1, -1, -1});
	local_ghosts[{}] += local_slice[{}] / 3.0;
}
void FillGhostCellsForLocalWithFineCornerNbr(const LocalData<3> &local_data, Corner<3> corner)
{
	auto local_slice  = local_data.getSliceOn(corner, {0, 0, 0});
	auto local_ghosts = local_data.getSliceOn(corner, {-1, -1, -1});
	local_ghosts[{}] += -local_slice[{}] / 3.0;
}
void FillLocalGhostCellsOnSides(const PatchInfo<3> &pinfo, const LocalData<3> &local_data)
{
	for (Side<3> side : Side<3>::getValues()) {
		if (pinfo.hasNbr(side)) {
			switch (pinfo.getNbrType(side)) {
				case NbrType::Normal:
					// do nothing
					break;
				case NbrType::Coarse:
					FillGhostCellsForLocalWithCoarseNbr(pinfo, local_data, side);
					break;
				case NbrType::Fine:
					FillGhostCellsForLocalWithFineNbr(local_data, side);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
void FillLocalGhostCellsOnEdges(const PatchInfo<3> &pinfo, const LocalData<3> &local_data)
{
	for (Edge edge : Edge::getValues()) {
		if (pinfo.hasNbr(edge)) {
			switch (pinfo.getNbrType(edge)) {
				case NbrType::Normal:
					// do nothing
					break;
				case NbrType::Coarse:
					FillGhostCellsForLocalWithCoarseEdgeNbr(pinfo, local_data, edge);
					break;
				case NbrType::Fine:
					FillGhostCellsForLocalWithFineEdgeNbr(local_data, edge);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
void FillLocalGhostCellsOnCorners(const PatchInfo<3> &pinfo, const LocalData<3> &local_data)
{
	for (Corner<3> corner : Corner<3>::getValues()) {
		if (pinfo.hasNbr(corner)) {
			switch (pinfo.getNbrType(corner)) {
				case NbrType::Normal:
					// do nothing
					break;
				case NbrType::Coarse:
					FillGhostCellsForLocalWithCoarseCornerNbr(pinfo, local_data, corner);
					break;
				case NbrType::Fine:
					FillGhostCellsForLocalWithFineCornerNbr(local_data, corner);
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
}
} // namespace

void TriLinearGhostFiller::fillGhostCellsForNbrPatch(const PatchInfo<3> &             pinfo,
                                                     const std::vector<LocalData<3>> &local_datas,
                                                     std::vector<LocalData<3>> &      nbr_datas,
                                                     Side<3>                          side,
                                                     NbrType                          nbr_type,
                                                     Orthant<2>                       orthant) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostCellsForNormalNbr(local_datas, nbr_datas, side);
			break;
		case NbrType::Coarse:
			FillGhostCellsForCoarseNbr(pinfo, local_datas, nbr_datas, side, orthant);
			break;
		case NbrType::Fine:
			FillGhostCellsForFineNbr(pinfo, local_datas, nbr_datas, side, orthant);
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void TriLinearGhostFiller::fillGhostCellsForEdgeNbrPatch(const PatchInfo<3> &             pinfo,
                                                         const std::vector<LocalData<3>> &local_datas,
                                                         std::vector<LocalData<3>> &      nbr_datas,
                                                         Edge                             edge,
                                                         NbrType                          nbr_type,
                                                         Orthant<1>                       orthant_on_coarse) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostCellsForNormalEdgeNbr(local_datas, nbr_datas, edge);
			break;
		case NbrType::Coarse:
			FillGhostCellsForCoarseEdgeNbr(local_datas, nbr_datas, edge, orthant_on_coarse);
			break;
		case NbrType::Fine:
			FillGhostCellsForFineEdgeNbr(local_datas, nbr_datas, edge, orthant_on_coarse);
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void TriLinearGhostFiller::fillGhostCellsForCornerNbrPatch(const PatchInfo<3> &             pinfo,
                                                           const std::vector<LocalData<3>> &local_datas,
                                                           std::vector<LocalData<3>> &      nbr_datas,
                                                           Corner<3>                        corner,
                                                           NbrType                          nbr_type) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostCellsForNormalCornerNbr(local_datas, nbr_datas, corner);
			break;
		case NbrType::Coarse:
			FillGhostCellsForCoarseCornerNbr(local_datas, nbr_datas, corner);
			break;
		case NbrType::Fine:
			FillGhostCellsForFineCornerNbr(local_datas, nbr_datas, corner);
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void TriLinearGhostFiller::fillGhostCellsForLocalPatch(const PatchInfo<3> &pinfo, std::vector<LocalData<3>> &local_datas) const
{
	for (const LocalData<3> &local_data : local_datas) {
		switch (this->getFillType()) {
			case GhostFillingType::Corners:
				FillLocalGhostCellsOnCorners(pinfo, local_data);
			case GhostFillingType::Edges:
				FillLocalGhostCellsOnEdges(pinfo, local_data);
			case GhostFillingType::Faces:
				FillLocalGhostCellsOnSides(pinfo, local_data);
				break;
			default:
				throw RuntimeError("Unsupported GhostFillingType");
		}
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