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
namespace
{
void FillGhostForNormalNbr(const std::vector<LocalData<2>> &local_datas,
                           const std::vector<LocalData<2>> &nbr_datas, const Side<2> side)
{
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto local_slice = local_datas[c].getSliceOnSide(side);
		auto nbr_ghosts  = nbr_datas[c].getGhostSliceOnSide(side.opposite(), 1);
		nested_loop<1>(
		nbr_ghosts.getStart(), nbr_ghosts.getEnd(),
		[&](const std::array<int, 1> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
	}
}
void FillGhostForCoarseNbr(std::shared_ptr<const PatchInfo<2>> pinfo,
                           const std::vector<LocalData<2>> &   local_datas,
                           const std::vector<LocalData<2>> &nbr_datas, const Side<2> side,
                           const Orthant<2> orthant)
{
	auto nbr_info = pinfo->getCoarseNbrInfo(side);
	int  offset   = 0;
	if (orthant.collapseOnAxis(side.getAxisIndex()) == Orthant<1>::upper()) {
		offset = pinfo->ns[!side.getAxisIndex()];
	}
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto local_slice = local_datas[c].getSliceOnSide(side);
		auto nbr_ghosts  = nbr_datas[c].getGhostSliceOnSide(side.opposite(), 1);
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               nbr_ghosts[{(coord[0] + offset) / 2}] += 2.0 / 3.0 * local_slice[coord];
		               });
	}
}
void FillGhostForFineNbr(std::shared_ptr<const PatchInfo<2>> pinfo,
                         const std::vector<LocalData<2>> &   local_datas,
                         const std::vector<LocalData<2>> &nbr_datas, const Side<2> side,
                         const Orthant<2> orthant)
{
	auto nbr_info = pinfo->getFineNbrInfo(side);
	int  offset   = 0;
	if (orthant.collapseOnAxis(side.getAxisIndex()) == Orthant<1>::upper()) {
		offset = pinfo->ns[!side.getAxisIndex()];
	}
	for (size_t c = 0; c < local_datas.size(); c++) {
		auto local_slice = local_datas[c].getSliceOnSide(side);
		auto nbr_ghosts  = nbr_datas[c].getGhostSliceOnSide(side.opposite(), 1);
		nested_loop<1>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               nbr_ghosts[coord] += 2.0 / 3.0 * local_slice[{(coord[0] + offset) / 2}];
		               });
	}
}
void FillLocalGhostsForCoarseNbr(std::shared_ptr<const PatchInfo<2>> pinfo,
                                 const LocalData<2> &local_data, const Side<2> side)
{
	auto local_slice  = local_data.getSliceOnSide(side);
	auto local_ghosts = local_data.getGhostSliceOnSide(side, 1);
	int  offset       = 0;
	if (pinfo->getCoarseNbrInfo(side).orth_on_coarse == Orthant<1>::upper()) {
		offset = pinfo->ns[!side.getAxisIndex()];
	}
	nested_loop<1>(local_ghosts.getStart(), local_ghosts.getEnd(),
	               [&](const std::array<int, 1> &coord) {
		               local_ghosts[coord] += 2.0 / 3.0 * local_slice[coord];
		               if ((coord[0] + offset) % 2 == 0) {
			               local_ghosts[{coord[0] + 1}] += -1.0 / 3.0 * local_slice[coord];
		               } else {
			               local_ghosts[{coord[0] - 1}] += -1.0 / 3.0 * local_slice[coord];
		               }
	               });
}
void FillLocalGhostsForFineNbr(const LocalData<2> &local_data, const Side<2> side)
{
	auto local_slice  = local_data.getSliceOnSide(side);
	auto local_ghosts = local_data.getGhostSliceOnSide(side, 1);
	nested_loop<1>(local_ghosts.getStart(), local_ghosts.getEnd(),
	               [&](const std::array<int, 1> &coord) {
		               local_ghosts[coord] += -1.0 / 3.0 * local_slice[coord];
	               });
}
} // namespace
void BiLinearGhostFiller::fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<2>> pinfo,
                                                    const std::vector<LocalData<2>> &   local_datas,
                                                    const std::vector<LocalData<2>> &   nbr_datas,
                                                    const Side<2> side, const NbrType nbr_type,
                                                    const Orthant<2> orthant) const
{
	switch (nbr_type) {
		case NbrType::Normal:
			FillGhostForNormalNbr(local_datas, nbr_datas, side);
			break;
		case NbrType::Coarse:
			FillGhostForCoarseNbr(pinfo, local_datas, nbr_datas, side, orthant);
			break;
		case NbrType::Fine:
			FillGhostForFineNbr(pinfo, local_datas, nbr_datas, side, orthant);
			break;
		default:
			throw RuntimeError("Unsupported Nbr Type");
	}
}

void BiLinearGhostFiller::fillGhostCellsForLocalPatch(
std::shared_ptr<const PatchInfo<2>> pinfo, const std::vector<LocalData<2>> &local_datas) const
{
	for (auto &local_data : local_datas) {
		for (Side<2> side : Side<2>::getValues()) {
			if (pinfo->hasNbr(side)) {
				switch (pinfo->getNbrType(side)) {
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
}
} // namespace ThunderEgg