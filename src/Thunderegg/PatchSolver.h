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

#ifndef THUNDEREGG_PATCHSOLVER_H
#define THUNDEREGG_PATCHSOLVER_H

#include <Thunderegg/Domain.h>
#include <Thunderegg/GMG/Smoother.h>
#include <Thunderegg/GhostFiller.h>
#include <Thunderegg/Operator.h>
#include <Thunderegg/Vector.h>

namespace Thunderegg
{
/**
 * @brief Solves the problem on the patches using a specified interface value
 *
 * @tparam D the number of cartesian dimensions
 */
template <size_t D> class PatchSolver : public virtual Operator<D>, public virtual GMG::Smoother<D>
{
	protected:
	std::shared_ptr<const Domain<D>>      domain;
	std::shared_ptr<const GhostFiller<D>> ghost_filler;

	public:
	/**
	 * @brief Destroy the Patch Solver object
	 */
	virtual ~PatchSolver() {}
	virtual void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                              const LocalData<D> f) const = 0;
	/**
	 * @brief Solve all the patches in the domain, assuming zero boundary conditions for the patches
	 *
	 * @param f the rhs vector
	 * @param u the lhs vector
	 */
	virtual void apply(std::shared_ptr<const Vector<D>> f,
	                   std::shared_ptr<Vector<D>>       u) const override
	{
		u->setWithGhost(0);
		for (std::shared_ptr<const PatchInfo<D>> pinfo : domain->getPatchInfoVector()) {
			solveSinglePatch(pinfo, u->getLocalData(pinfo->local_index),
			                 f->getLocalData(pinfo->local_index));
		}
	}
	/**
	 * @brief Solve all the patches in the domain, using the values in u for the boundary conditions
	 *
	 * @param f the rhs vector
	 * @param u the lhs vector
	 */
	virtual void smooth(std::shared_ptr<const Vector<D>> f,
	                    std::shared_ptr<Vector<D>>       u) const override
	{
		ghost_filler->fillGhost(u);
		for (std::shared_ptr<const PatchInfo<D>> pinfo : domain->getPatchInfoVector()) {
			solveSinglePatch(pinfo, u->getLocalData(pinfo->local_index),
			                 f->getLocalData(pinfo->local_index));
		}
	}
};
} // namespace Thunderegg
#endif
