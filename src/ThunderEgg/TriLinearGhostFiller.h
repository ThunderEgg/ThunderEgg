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
 */
class TriLinearGhostFiller : public MPIGhostFiller<3>
{
	private:
	void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<3>> pinfo,
	                               const LocalData<3> local_data, const LocalData<3> nbr_data,
	                               const Side<3> side, const NbrType nbr_type,
	                               const Orthant<3> orthant) const override;

	void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<3>> pinfo,
	                                 const LocalData<3>                  local_data) const override;

	public:
	/**
	 * @brief Construct a new TriLinearGhostFiller object
	 *
	 * @param domain the domain on which ghosts will be filled
	 */
	TriLinearGhostFiller(std::shared_ptr<const Domain<3>> domain) : MPIGhostFiller<3>(domain, 1) {}
};
} // namespace ThunderEgg
#endif
