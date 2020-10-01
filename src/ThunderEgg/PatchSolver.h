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

#ifndef THUNDEREGG_PATCHSOLVER_H
#define THUNDEREGG_PATCHSOLVER_H

#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/Smoother.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/Operator.h>
#include <ThunderEgg/Timer.h>
#include <ThunderEgg/Vector.h>

namespace ThunderEgg
{
/**
 * @brief Solves the problem on the patches using a specified interface value
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D> class PatchSolver : public virtual Operator<D>, public virtual GMG::Smoother<D>
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
	 * @brief Construct a new PatchSolver object
	 *
	 * 	This sets the Domain and the GhostFiller that the object uses.
	 *
	 * @param domain the Domain
	 * @param ghost_filler the GhostFiller
	 */
	PatchSolver(std::shared_ptr<const Domain<D>>      domain,
	            std::shared_ptr<const GhostFiller<D>> ghost_filler)
	: domain(domain), ghost_filler(ghost_filler)
	{
	}
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
	 * @brief Get the Domain object
	 *
	 * @return std::shared_ptr<const Domain<D>> the Domain
	 */
	std::shared_ptr<const Domain<D>> getDomain() const
	{
		return domain;
	}
	/**
	 * @brief Get the GhostFiller object
	 *
	 * @return std::shared_ptr<const GhostFiller<D>> the GhostFiller
	 */
	std::shared_ptr<const GhostFiller<D>> getGhostFiller() const
	{
		return ghost_filler;
	}
	/**
	 * @brief Perform a single solve over a patch
	 *
	 * @param pinfo the PatchInfo for the patch
	 * @param us the left hand side
	 * @param fs the right hand side
	 */
	virtual void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                              const std::vector<LocalData<D>> &   fs,
	                              std::vector<LocalData<D>> &         us) const = 0;
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
			auto us = u->getLocalDatas(pinfo->local_index);
			auto fs = f->getLocalDatas(pinfo->local_index);
			solveSinglePatch(pinfo, us, fs);
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
			auto us = u->getLocalDatas(pinfo->local_index);
			auto fs = f->getLocalDatas(pinfo->local_index);
			solveSinglePatch(pinfo, us, fs);
			if (timer) {
				timer->stop("Single Patch Solve");
			}
		}
		if (timer) {
			timer->stop("Total Patch Smooth");
		}
	}
};
} // namespace ThunderEgg
#endif
