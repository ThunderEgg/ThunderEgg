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

#ifndef THUNDEREGG_PATCHOPERATOR_H
#define THUNDEREGG_PATCHOPERATOR_H
#include <Thunderegg/GhostFiller.h>
#include <Thunderegg/Operator.h>
#include <Thunderegg/Vector.h>
namespace Thunderegg
{
template <size_t D> class PatchOperator : public Operator<D>
{
	protected:
	/**
	 * @brief the domain that is being solved over
	 */
	std::shared_ptr<const Domain<D>> domain;
	/**
	 * @brief The ghost filler, needed for smoothing
	 */
	std::shared_ptr<const GhostFiller<D>> ghost_filler;

	public:
	virtual ~PatchOperator() {}

	virtual void applySinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                              LocalData<D> f) const = 0;
	virtual void addGhostToRHS(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                           LocalData<D> f) const    = 0;

	virtual void apply(std::shared_ptr<const Vector<D>> u,
	                   std::shared_ptr<Vector<D>>       f) const override
	{
		ghost_filler->fillGhost(u);
		for (auto pinfo : domain->getPatchInfoVector()) 
        {
			applySinglePatch(pinfo, u->getLocalData(pinfo->local_index),
			                 f->getLocalData(pinfo->local_index));
		}
	}
};
} // namespace Thunderegg
#endif
