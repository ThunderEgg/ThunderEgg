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

#ifndef THUNDEREGG_SIMPLEGHOSTFILLER_H
#define THUNDEREGG_SIMPLEGHOSTFILLER_H
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/Vector.h>
namespace ThunderEgg
{
/**
 * @brief Exchanges ghost cells on patches
 *
 * @tparam D the number of Cartesian dimensions in the patches.
 */
template <int D> class SimpleGhostFiller : public GhostFiller<D>
{
	protected:
	std::shared_ptr<const Domain<D>> domain;

	public:
	SimpleGhostFiller(std::shared_ptr<const Domain<D>> domain)
	{
		this->domain = domain;
	}

	/**
	 * @brief Fill the ghost cells on a vector
	 *
	 * @param u  the vector
	 */
	void fillGhost(std::shared_ptr<const Vector<D>> u) const
	{
		for (auto pinfo : domain->getPatchInfoVector()) {
			for (Side<D> s : Side<D>::getValues()) {
				if (pinfo->hasNbr(s)) {
					auto               nbr_info    = pinfo->getNormalNbrInfo(s);
					const LocalData<D> this_patch  = u->getLocalData(pinfo->local_index);
					const LocalData<D> other_patch = u->getLocalData(nbr_info.local_index);
					auto               this_side   = this_patch.getSliceOnSide(s);
					auto other_side = other_patch.getGhostSliceOnSide(s.opposite(), 1);
					nested_loop<D - 1>(this_side.getStart(), this_side.getEnd(),
					                   [&](const std::array<int, D - 1> &coord) {
						                   other_side[coord] = this_side[coord];
					                   });
				}
			}
		}
	}
};
} // namespace ThunderEgg
#endif
