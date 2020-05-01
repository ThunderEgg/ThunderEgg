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
#include <Thunderegg/Timer.h>
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
	/**
	 * @brief the domain that is being solved over
	 */
	std::shared_ptr<const Domain<D>> domain;
	/**
	 * @brief The ghost filler, needed for smoothing
	 */
	std::shared_ptr<const GhostFiller<D>> ghost_filler;
	/**
	 * @brief The timer
	 */
	mutable std::shared_ptr<Timer> timer;

	public:
	/**
	 * @brief Destroy the Patch Solver object
	 */
	virtual ~PatchSolver() {}
	/**
	 * @brief Set the Timer object
	 *
	 * @param timer the timer
	 */
	void setTimer(std::shared_ptr<Timer> timer) const
	{
		this->timer = timer;
	}
	/**
	 * @brief Get the Timer object
	 *
	 * @return std::shared_ptr<Timer> the timer
	 */
	std::shared_ptr<Timer> getTimer() const
	{
		return timer;
	}
	/**
	 * @brief Perform a single solve over a patch
	 *
	 * @param pinfo the PatchInfo for the patch
	 * @param u the left hand side
	 * @param f the right hand side
	 */
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
		if (timer) {
			timer->start("Total Patch Solve");
		}
		for (std::shared_ptr<const PatchInfo<D>> pinfo : domain->getPatchInfoVector()) {
			if (timer) {
				timer->start("Single Patch Solve");
			}
			solveSinglePatch(pinfo, u->getLocalData(pinfo->local_index),
			                 f->getLocalData(pinfo->local_index));
			if (timer) {
				timer->stop("Single Patch Solve");
			}
		}
		if (timer) {
			timer->stop("Total Patch Solve");
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
		if (timer) {
			timer->start("Total Patch Smooth");
		}
		ghost_filler->fillGhost(u);
		for (std::shared_ptr<const PatchInfo<D>> pinfo : domain->getPatchInfoVector()) {
			if (timer) {
				timer->start("Single Patch Solve");
			}
			solveSinglePatch(pinfo, u->getLocalData(pinfo->local_index),
			                 f->getLocalData(pinfo->local_index));
			if (timer) {
				timer->stop("Single Patch Solve");
			}
		}
		if (timer) {
			timer->stop("Total Patch Smooth");
		}
	}
};
} // namespace Thunderegg
#endif
