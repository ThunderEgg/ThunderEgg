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
/**
 * @brief Get the Offset on a given orthant
 *
 * @param ns the dimensions of the patch
 * @param orth the orthant
 * @return std::array<int,2> the offsets
 */
static std::array<int, 2> getOffset(const std::array<int, 3> ns, Side<3> s, Orthant<2> orth)
{
	std::array<int, 2> offset = {0, 0};
	for (size_t i = 0; i < s.getAxisIndex(); i++) {
		if (orth.isOnSide(Side<2>::HigherSideOnAxis(i))) {
			offset[i] = ns[i];
		}
	}
	for (size_t i = s.getAxisIndex() + 1; i < 3; i++) {
		if (orth.isOnSide(Side<2>::HigherSideOnAxis(i - 1))) {
			offset[i - 1] = ns[i];
		}
	}
	return offset;
}
void TriLinearGhostFiller::fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<3>> pinfo,
                                                     const LocalData<3> &                local_data,
                                                     const LocalData<3> &                nbr_data,
                                                     const Side<3> side, const NbrType nbr_type,
                                                     const Orthant<3> orthant) const
{
	if (nbr_type == NbrType::Normal) {
		auto local_slice = local_data.getSliceOnSide(side);
		auto nbr_ghosts  = nbr_data.getGhostSliceOnSide(side.opposite(), 1);
		nested_loop<2>(
		nbr_ghosts.getStart(), nbr_ghosts.getEnd(),
		[&](const std::array<int, 2> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
	} else if (nbr_type == NbrType::Coarse) {
		auto               nbr_info    = pinfo->getCoarseNbrInfo(side);
		auto               local_slice = local_data.getSliceOnSide(side);
		auto               nbr_ghosts  = nbr_data.getGhostSliceOnSide(side.opposite(), 1);
		std::array<int, 2> offset
		= getOffset(pinfo->ns, side, orthant.collapseOnAxis(side.getAxisIndex()));
		nested_loop<2>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(),
		               [&](const std::array<int, 2> &coord) {
			               std::array<int, 2> coarse_coord;
			               for (int i = 0; i < 2; i++) {
				               coarse_coord[i] = (coord[i] + offset[i]) / 2;
			               }
			               nbr_ghosts[coarse_coord] += 1.0 / 3.0 * local_slice[coord];
		               });
	} else if (nbr_type == NbrType::Fine) {
		auto               nbr_info    = pinfo->getFineNbrInfo(side);
		auto               local_slice = local_data.getSliceOnSide(side);
		auto               nbr_ghosts  = nbr_data.getGhostSliceOnSide(side.opposite(), 1);
		std::array<int, 2> offset
		= getOffset(pinfo->ns, side, orthant.collapseOnAxis(side.getAxisIndex()));
		nested_loop<2>(nbr_ghosts.getStart(), nbr_ghosts.getEnd(),
		               [&](const std::array<int, 2> &coord) {
			               std::array<int, 2> coarse_coord;
			               for (int i = 0; i < 2; i++) {
				               coarse_coord[i] = (coord[i] + offset[i]) / 2;
			               }
			               nbr_ghosts[coord] += 4.0 / 6.0 * local_slice[coarse_coord];
		               });
	}
}

void TriLinearGhostFiller::fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<3>> pinfo,
                                                       const LocalData<3> &local_data) const
{
	for (Side<3> side : Side<3>::getValues()) {
		if (pinfo->hasNbr(side)) {
			NbrType nbr_type = pinfo->getNbrType(side);
			if (nbr_type == NbrType::Coarse) {
				auto               local_slice  = local_data.getSliceOnSide(side);
				auto               local_ghosts = local_data.getGhostSliceOnSide(side, 1);
				auto               nbr_info     = pinfo->getCoarseNbrInfo(side);
				std::array<int, 2> offset = getOffset(pinfo->ns, side, nbr_info.orth_on_coarse);
				nested_loop<2>(
				local_ghosts.getStart(), local_ghosts.getEnd(),
				[&](const std::array<int, 2> &coord) {
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
			} else if (nbr_type == NbrType::Fine) {
				auto local_slice  = local_data.getSliceOnSide(side);
				auto local_ghosts = local_data.getGhostSliceOnSide(side, 1);
				nested_loop<2>(local_ghosts.getStart(), local_ghosts.getEnd(),
				               [&](const std::array<int, 2> &coord) {
					               local_ghosts[coord] -= 1.0 / 3.0 * local_slice[coord];
				               });
			}
		}
	}
}
TriLinearGhostFiller::TriLinearGhostFiller(std::shared_ptr<const Domain<3>> domain)
: MPIGhostFiller<3>(domain, 1)
{
	for (int n : domain->getNs()) {
		if (n % 2 != 0) {
			throw RuntimeError("TriLinearGhostFiller only supports even number patch sizes");
		}
	}
}

} // namespace ThunderEgg