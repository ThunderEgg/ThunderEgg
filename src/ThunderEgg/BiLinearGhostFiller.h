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

#ifndef THUNDEREGG_BILINEARGHOSTFILLER_H
#define THUNDEREGG_BILINEARGHOSTFILLER_H
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/MPIGhostFiller.h>
namespace ThunderEgg
{
/**
 * @brief Exchanges ghost cells on patches, uses a BiLinear interpolation scheme for refinement
 * boundaries
 */
class BiLinearGhostFiller : public MPIGhostFiller<2>
{
	public:
	/**
	 * @brief Construct a new BiLinearGhostFiller object
	 *
	 * @param domain_in the domain to fill ghosts for
	 */
	BiLinearGhostFiller(std::shared_ptr<const Domain<2>> domain_in) : MPIGhostFiller<2>(domain_in, 1) {}
	void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<2>> pinfo,
	                               const std::vector<LocalData<2>> &   local_datas,
	                               const std::vector<LocalData<2>> &   nbr_datas,
	                               Side<2>                             sides,
	                               NbrType                             nbr_type,
	                               Orthant<1>                          orthant_on_coarse) const override;

	void fillGhostCellsForEdgeNbrPatch(std::shared_ptr<const PatchInfo<2>> pinfo,
	                                   const std::vector<LocalData<2>> &   local_datas,
	                                   const std::vector<LocalData<2>> &   nbr_datas,
	                                   Edge<2>                             edge,
	                                   NbrType                             nbr_type,
	                                   Orthant<1>                          orthant_on_coarse) const override
	{
	}

	void fillGhostCellsForCornerNbrPatch(std::shared_ptr<const PatchInfo<2>> pinfo,
	                                     const std::vector<LocalData<2>> &   local_datas,
	                                     const std::vector<LocalData<2>> &   nbr_datas,
	                                     Corner<2>                           corner,
	                                     NbrType                             nbr_type) const override
	{
	}
	void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<2>> pinfo, const std::vector<LocalData<2>> &local_datas) const override;
};
} // namespace ThunderEgg
#endif
