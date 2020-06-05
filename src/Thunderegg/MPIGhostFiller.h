/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_MPIGHOSTFILLER_H
#define THUNDEREGG_MPIGHOSTFILLER_H

#include <Thunderegg/GhostFiller.h>
namespace Thunderegg
{
template <size_t D> class MPIGhostFiller : public GhostFiller<D>
{
	private:
	virtual void fillGhostCellsForNbrPatch(const LocalData<D> local_data, LocalData<D> nbr_data,
	                                       const std::vector<Side<D>> &side, const NbrType nbr_type,
	                                       const Orthant<D> orthant) const = 0;

	virtual void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                         LocalData<D> local_data) const = 0;

	protected:
	std::shared_ptr<Domain<D>> domain;
	int                        side_cases;

	public:
	MPIGhostFiller(std::shared_ptr<Domain<D>> domain_in, int side_cases_in)
	: domain(domain_in), side_cases(side_cases_in)
	{
	}
	/**
	 * @brief Fill ghost cells on a vector
	 *
	 * @param u  the vector
	 */
	void fillGhost(std::shared_ptr<const Vector<D>> u) const {}
};
} // namespace Thunderegg
#endif