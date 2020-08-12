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

#ifndef THUNDEREGG_BIQUADRATICGHOSTFILLER_H
#define THUNDEREGG_BIQUADRATICGHOSTFILLER_H
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/MPIGhostFiller.h>
namespace ThunderEgg
{
/**
 * @brief Exchanges ghost cells on patches
 */
class BiQuadraticGhostFiller : public MPIGhostFiller<2>
{
	private:
	void fillGhostForLocalWithCoarseNbr(const LocalData<1> &inner_slice, const LocalData<1> &slice,
	                                    LocalData<1> &ghost) const
	{
		int n = ghost.getLengths()[0];
		for (int idx = 0; idx < n; idx++) {
			ghost[{idx}] += 2 * slice[{idx}] / 3 - inner_slice[{idx}] / 5;
		}
	}
	void fillGhostForLocalWithFineNbr(const LocalData<1> &slice, LocalData<1> &ghost) const
	{
		int n = ghost.getLengths()[0];
		ghost[{0}] += -slice[{0}] / 10 + slice[{1}] / 15 - slice[{2}] / 30;
		for (int idx = 1; idx < n - 1; idx++) {
			ghost[{idx}] += -slice[{idx - 1}] / 30 - slice[{idx + 1}] / 30;
		}
		ghost[{n - 1}] += -slice[{n - 1}] / 10 + slice[{n - 2}] / 15 - slice[{n - 3}] / 30;
	}
	void fillGhostForCoarseNbrLower(const LocalData<1> &inner_slice, const LocalData<1> &slice,
	                                LocalData<1> &ghost) const
	{
		int n = ghost.getLengths()[0];
		for (int idx = 0; idx < n; idx++) {
			ghost[{idx / 2}] += slice[{idx}] / 3 + inner_slice[{idx}] / 5;
		}
	}
	void fillGhostForCoarseNbrUpper(const LocalData<1> &inner_slice, const LocalData<1> &slice,
	                                LocalData<1> &ghost) const
	{
		int n = ghost.getLengths()[0];
		for (int idx = 0; idx < n; idx++) {
			ghost[{(idx + n) / 2}] += slice[{idx}] / 3 + inner_slice[{idx}] / 5;
		}
	}
	void fillGhostForFineNbrLower(const LocalData<1> &slice, LocalData<1> &ghost) const
	{
		int n = ghost.getLengths()[0];
		ghost[{0}] += 3 * slice[{0}] / 4 - 3 * slice[{1}] / 10 + slice[{2}] / 12;
		ghost[{1}] += 7 * slice[{0}] / 20 + 7 * slice[{1}] / 30 - slice[{2}] / 20;
		for (int idx = 2; idx < n; idx++) {
			if (idx % 2 == 0) {
				ghost[{idx}]
				+= slice[{idx / 2 - 1}] / 12 + slice[{idx / 2}] / 2 - slice[{idx / 2 + 1}] / 20;
			} else {
				ghost[{idx}]
				+= -slice[{idx / 2 - 1}] / 20 + slice[{idx / 2}] / 2 + slice[{idx / 2 + 1}] / 12;
			}
		}
	}
	void fillGhostForFineNbrUpper(const LocalData<1> &slice, LocalData<1> &ghost) const
	{
		int n = ghost.getLengths()[0];
		for (int idx = 0; idx < n - 2; idx++) {
			if ((idx + n) % 2 == 0) {
				ghost[{idx}] += slice[{(idx + n) / 2 - 1}] / 12 + slice[{(idx + n) / 2}] / 2
				                - slice[{(idx + n) / 2 + 1}] / 20;
			} else {
				ghost[{idx}] += -slice[{(idx + n) / 2 - 1}] / 20 + slice[{(idx + n) / 2}] / 2
				                + slice[{(idx + n) / 2 + 1}] / 12;
			}
		}
		ghost[{n - 2}] += 7 * slice[{n - 1}] / 20 + 7 * slice[{n - 2}] / 30 - slice[{n - 3}] / 20;
		ghost[{n - 1}] += 3 * slice[{n - 1}] / 4 - 3 * slice[{n - 2}] / 10 + slice[{n - 3}] / 12;
	}
	void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<2>> pinfo,
	                               const LocalData<2> local_data, const LocalData<2> nbr_data,
	                               const Side<2> side, const NbrType nbr_type,
	                               const Orthant<2> orthant) const override
	{
		switch (nbr_type) {
			case NbrType::Normal: {
				auto local_slice = local_data.getSliceOnSide(side);
				auto nbr_ghosts  = nbr_data.getGhostSliceOnSide(side.opposite(), 1);
				nested_loop<1>(
				nbr_ghosts.getStart(), nbr_ghosts.getEnd(),
				[&](const std::array<int, 1> &coord) { nbr_ghosts[coord] = local_slice[coord]; });
			} break;
			case NbrType::Coarse: {
				auto local_slice       = local_data.getSliceOnSide(side);
				auto inner_local_slice = local_data.getSliceOnSide(side, 1);
				auto nbr_ghosts        = nbr_data.getGhostSliceOnSide(side.opposite(), 1);
				if (orthant.collapseOnAxis(side.getAxisIndex()) == Orthant<1>::lower()) {
					fillGhostForCoarseNbrLower(inner_local_slice, local_slice, nbr_ghosts);
				} else {
					fillGhostForCoarseNbrUpper(inner_local_slice, local_slice, nbr_ghosts);
				}

			} break;
			case NbrType::Fine: {
				auto local_slice = local_data.getSliceOnSide(side);
				auto nbr_ghosts  = nbr_data.getGhostSliceOnSide(side.opposite(), 1);
				if (orthant.collapseOnAxis(side.getAxisIndex()) == Orthant<1>::lower()) {
					fillGhostForFineNbrLower(local_slice, nbr_ghosts);
				} else {
					fillGhostForFineNbrUpper(local_slice, nbr_ghosts);
				}
			} break;
		}
	}

	void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<2>> pinfo,
	                                 const LocalData<2>                  local_data) const override
	{
		for (Side<2> side : Side<2>::getValues()) {
			if (pinfo->hasNbr(side)) {
				switch (pinfo->getNbrType(side)) {
					case NbrType::Normal:
						// nothing need to be done
						break;
					case NbrType::Coarse: {
						auto inner_local_slice = local_data.getSliceOnSide(side, 1);
						auto local_slice       = local_data.getSliceOnSide(side);
						auto local_ghosts      = local_data.getGhostSliceOnSide(side, 1);
						fillGhostForLocalWithCoarseNbr(inner_local_slice, local_slice,
						                               local_ghosts);
					} break;
					case NbrType::Fine: {
						auto local_slice  = local_data.getSliceOnSide(side);
						auto local_ghosts = local_data.getGhostSliceOnSide(side, 1);
						fillGhostForLocalWithFineNbr(local_slice, local_ghosts);
					} break;
				}
			}
		}
	}

	public:
	BiQuadraticGhostFiller(std::shared_ptr<const Domain<2>> domain_in)
	: MPIGhostFiller<2>(domain_in, 1)
	{
	}
};
} // namespace ThunderEgg
#endif
