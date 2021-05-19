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
 * @brief Exchanges ghost cells on patches, handles refinement boundaries with a biquadratic
 * interpolation scheme
 */
class BiQuadraticGhostFiller : public MPIGhostFiller<2>
{
	public:
	/**
	 * @brief Construct a new BiQuadraticGhostFiller object
	 *
	 * @param domain the domain that is being fill for
	 * @param fill_type the ghost filling type to perform
	 */
	BiQuadraticGhostFiller(std::shared_ptr<const Domain<2>> domain, GhostFillingType fill_type);

	void fillGhostCellsForNbrPatch(const PatchInfo<2> &             pinfo,
	                               const std::vector<LocalData<2>> &local_datas,
	                               std::vector<LocalData<2>> &      nbr_datas,
	                               Side<2>                          side,
	                               NbrType                          nbr_type,
	                               Orthant<1>                       orthant_on_coarse) const override;

	void fillGhostCellsForEdgeNbrPatch(const PatchInfo<2> &             pinfo,
	                                   const std::vector<LocalData<2>> &local_datas,
	                                   std::vector<LocalData<2>> &      nbr_datas,
	                                   Edge                             edge,
	                                   NbrType                          nbr_type,
	                                   Orthant<1>                       orthant_on_coarse) const override;

	void fillGhostCellsForCornerNbrPatch(const PatchInfo<2> &             pinfo,
	                                     const std::vector<LocalData<2>> &local_datas,
	                                     std::vector<LocalData<2>> &      nbr_datas,
	                                     Corner<2>                        corner,
	                                     NbrType                          nbr_type) const override;

	void fillGhostCellsForLocalPatch(const PatchInfo<2> &pinfo, std::vector<LocalData<2>> &local_datas) const override;
};
} // namespace ThunderEgg
#endif
