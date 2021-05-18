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

#include <ThunderEgg/BiQuadraticGhostFiller.h>
#include <ThunderEgg/RuntimeError.h>
namespace ThunderEgg
{
BiQuadraticGhostFiller::BiQuadraticGhostFiller(std::shared_ptr<const Domain<2>> domain_in) : MPIGhostFiller<2>(domain_in, GhostFillingType::Faces) {}

namespace
{
void FillGhostForLocalWithCoarseNbr(const LocalData<2> &local_data, const Side<2> side)
{
	auto inner_slice = local_data.getSliceOn(side, {1});
	auto slice       = local_data.getSliceOn(side, {0});
	auto ghost       = local_data.getSliceOn(side, {-1});
	int  n           = ghost.getLengths()[0];
	for (int idx = 0; idx < n; idx++) {
		ghost[{idx}] += 2 * slice[{idx}] / 3 - inner_slice[{idx}] / 5;
	}
}
void FillGhostForLocalWithFineNbr(const LocalData<2> &local_data, const Side<2> side)
{
	auto slice = local_data.getSliceOn(side, {0});
	auto ghost = local_data.getSliceOn(side, {-1});
	int  n     = ghost.getLengths()[0];
	ghost[{0}] += -slice[{0}] / 10 + slice[{1}] / 15 - slice[{2}] / 30;
	for (int idx = 1; idx < n - 1; idx++) {
		ghost[{idx}] += -slice[{idx - 1}] / 30 - slice[{idx + 1}] / 30;
	}
	ghost[{n - 1}] += -slice[{n - 1}] / 10 + slice[{n - 2}] / 15 - slice[{n - 3}] / 30;
}
void FillGhostForNormalNbr(const std::vector<LocalData<2>> &local_datas, std::vector<LocalData<2>> &nbr_datas, const Side<2> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto local_slice = local_datas[c].getSliceOn(side, {0});
		auto nbr_ghosts  = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(), [&](const std::array<int, 1> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
	}
}
void FillGhostForCoarseNbrLower(const std::vector<LocalData<2>> &local_datas, std::vector<LocalData<2>> &nbr_datas, const Side<2> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto slice       = local_datas[c].getSliceOn(side, {0});
		auto inner_slice = local_datas[c].getSliceOn(side, {1});
		auto ghost       = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		int  n           = ghost.getLengths()[0];
		for (int idx = 0; idx < n; idx++) {
			ghost[{idx / 2}] += slice[{idx}] / 3 + inner_slice[{idx}] / 5;
		}
	}
}
void FillGhostForCoarseNbrUpper(const std::vector<LocalData<2>> &local_datas, std::vector<LocalData<2>> &nbr_datas, const Side<2> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto slice       = local_datas[c].getSliceOn(side, {0});
		auto inner_slice = local_datas[c].getSliceOn(side, {1});
		auto ghost       = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		int  n           = ghost.getLengths()[0];
		for (int idx = 0; idx < n; idx++) {
			ghost[{(idx + n) / 2}] += slice[{idx}] / 3 + inner_slice[{idx}] / 5;
		}
	}
}
void FillGhostForFineNbrLower(const std::vector<LocalData<2>> &local_datas, std::vector<LocalData<2>> &nbr_datas, const Side<2> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto slice = local_datas[c].getSliceOn(side, {0});
		auto ghost = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		int  n     = ghost.getLengths()[0];
		ghost[{0}] += 3 * slice[{0}] / 4 - 3 * slice[{1}] / 10 + slice[{2}] / 12;
		ghost[{1}] += 7 * slice[{0}] / 20 + 7 * slice[{1}] / 30 - slice[{2}] / 20;
		for (int idx = 2; idx < n; idx++) {
			if (idx % 2 == 0) {
				ghost[{idx}] += slice[{idx / 2 - 1}] / 12 + slice[{idx / 2}] / 2 - slice[{idx / 2 + 1}] / 20;
			} else {
				ghost[{idx}] += -slice[{idx / 2 - 1}] / 20 + slice[{idx / 2}] / 2 + slice[{idx / 2 + 1}] / 12;
			}
		}
	}
}
void FillGhostForFineNbrUpper(const std::vector<LocalData<2>> &local_datas, std::vector<LocalData<2>> &nbr_datas, const Side<2> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto slice = local_datas[c].getSliceOn(side, {0});
		auto ghost = nbr_datas[c].getSliceOn(side.opposite(), {-1});
		int  n     = ghost.getLengths()[0];
		for (int idx = 0; idx < n - 2; idx++) {
			if ((idx + n) % 2 == 0) {
				ghost[{idx}] += slice[{(idx + n) / 2 - 1}] / 12 + slice[{(idx + n) / 2}] / 2 - slice[{(idx + n) / 2 + 1}] / 20;
			} else {
				ghost[{idx}] += -slice[{(idx + n) / 2 - 1}] / 20 + slice[{(idx + n) / 2}] / 2 + slice[{(idx + n) / 2 + 1}] / 12;
			}
		}
		ghost[{n - 2}] += 7 * slice[{n - 1}] / 20 + 7 * slice[{n - 2}] / 30 - slice[{n - 3}] / 20;
		ghost[{n - 1}] += 3 * slice[{n - 1}] / 4 - 3 * slice[{n - 2}] / 10 + slice[{n - 3}] / 12;
	}
}
} // namespace

void BiQuadraticGhostFiller::fillGhostCellsForNbrPatch(const PatchInfo<2> &             pinfo,
                                                       const std::vector<LocalData<2>> &local_datas,
                                                       std::vector<LocalData<2>> &      nbr_datas,
                                                       Side<2>                          side,
                                                       NbrType                          nbr_type,
                                                       Orthant<1>                       orthant) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostForNormalNbr(local_datas, nbr_datas, side);
			break;
		case NbrType::Coarse:
			if (orthant == Orthant<1>::lower()) {
				FillGhostForCoarseNbrLower(local_datas, nbr_datas, side);
			} else {
				FillGhostForCoarseNbrUpper(local_datas, nbr_datas, side);
			}
			break;
		case NbrType::Fine:
			if (orthant == Orthant<1>::lower()) {
				FillGhostForFineNbrLower(local_datas, nbr_datas, side);
			} else {
				FillGhostForFineNbrUpper(local_datas, nbr_datas, side);
			}
			break;
		default:
			throw RuntimeError("Unsupported NbrType");
	}
}

void BiQuadraticGhostFiller::fillGhostCellsForEdgeNbrPatch(const PatchInfo<2> &             pinfo,
                                                           const std::vector<LocalData<2>> &local_datas,
                                                           std::vector<LocalData<2>> &      nbr_datas,
                                                           Edge                             edge,
                                                           NbrType                          nbr_type,
                                                           Orthant<1>                       orthant_on_coarse) const
{
	// no edges for 2d
}

void BiQuadraticGhostFiller::fillGhostCellsForCornerNbrPatch(const PatchInfo<2> &             pinfo,
                                                             const std::vector<LocalData<2>> &local_datas,
                                                             std::vector<LocalData<2>> &      nbr_datas,
                                                             Corner<2>                        corner,
                                                             NbrType                          nbr_type) const
{
	// corners not implimented
}

void BiQuadraticGhostFiller::fillGhostCellsForLocalPatch(const PatchInfo<2> &pinfo, std::vector<LocalData<2>> &local_datas) const
{
	for (const LocalData<2> &local_data : local_datas) {
		for (Side<2> side : Side<2>::getValues()) {
			if (pinfo.hasNbr(side)) {
				switch (pinfo.getNbrType(side)) {
					case NbrType::Normal:
						// nothing need to be done
						break;
					case NbrType::Coarse:
						FillGhostForLocalWithCoarseNbr(local_data, side);
						break;
					case NbrType::Fine:
						FillGhostForLocalWithFineNbr(local_data, side);
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
		}
	}
}
} // namespace ThunderEgg
