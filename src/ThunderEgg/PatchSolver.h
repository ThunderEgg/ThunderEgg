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
	private:
	/**
	 * @brief the domain that is being solved over
	 */
	Domain<D> domain;
	/**
	 * @brief The ghost filler, needed for smoothing
	 */
	std::shared_ptr<const GhostFiller<D>> ghost_filler;

	public:
	/**
	 * @brief Construct a new PatchSolver object
	 *
	 * 	This sets the Domain and the GhostFiller that the object uses.
	 *
	 * @param domain the Domain
	 * @param ghost_filler the GhostFiller
	 */
	PatchSolver(const Domain<D> &domain, const GhostFiller<D> &ghost_filler) : domain(domain), ghost_filler(ghost_filler.clone()) {}
	/**
	 * @brief Destroy the Patch Solver object
	 */
	virtual ~PatchSolver() {}
	/**
	 * @brief Clone this patch solver
	 *
	 * @return PatchSolver<D>* a newly allocated copy of this patch solver
	 */
	virtual PatchSolver<D> *clone() const override = 0;
	/**
	 * @brief Get the Domain object
	 *
	 * @return const Domain<D>& the Domain
	 */
	const Domain<D> &getDomain() const
	{
		return domain;
	}
	/**
	 * @brief Get the GhostFiller object
	 *
	 * @return const GhostFiller<D>& the GhostFiller
	 */
	const GhostFiller<D> &getGhostFiller() const
	{
		return *ghost_filler;
	}
	/**
	 * @brief Perform a single solve over a patch
	 *
	 * @param pinfo the PatchInfo for the patch
	 * @param f_view the left hand side
	 * @param u_view the right hand side
	 */
	virtual void solveSinglePatch(const PatchInfo<D> &pinfo, const PatchView<const double, D> &f_view, const PatchView<double, D> &u_view) const = 0;
	/**
	 * @brief Solve all the patches in the domain, assuming zero boundary conditions for the patches
	 *
	 * @param f the rhs vector
	 * @param u the lhs vector
	 */
	virtual void apply(const Vector<D> &f, Vector<D> &u) const override
	{
		if constexpr (ENABLE_DEBUG) {
			if (u.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
				throw RuntimeError("u vector is incorrect length");
			}
			if (f.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
				throw RuntimeError("f vector is incorrect length");
			}
		}
		u.setWithGhost(0);
		if (domain.hasTimer()) {
			domain.getTimer()->startDomainTiming(domain.getId(), "Total Patch Solve");
		}
		for (const PatchInfo<D> &pinfo : domain.getPatchInfoVector()) {
			if (domain.hasTimer()) {
				domain.getTimer()->start("Single Patch Solve");
			}
			PatchView<const double, D> f_view = f.getPatchView(pinfo.local_index);
			PatchView<double, D>       u_view = u.getPatchView(pinfo.local_index);
			solveSinglePatch(pinfo, f_view, u_view);
			if (domain.hasTimer()) {
				domain.getTimer()->stop("Single Patch Solve");
			}
		}
		if (domain.hasTimer()) {
			domain.getTimer()->stopDomainTiming(domain.getId(), "Total Patch Solve");
		}
	}
	/**
	 * @brief Solve all the patches in the domain, using the values in u for the boundary conditions
	 *
	 * @param f the rhs vector
	 * @param u the lhs vector
	 */
	virtual void smooth(const Vector<D> &f, Vector<D> &u) const override
	{
		if constexpr (ENABLE_DEBUG) {
			if (u.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
				throw RuntimeError("u vector is incorrect length");
			}
			if (f.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
				throw RuntimeError("f vector is incorrect length");
			}
		}
		if (domain.hasTimer()) {
			domain.getTimer()->startDomainTiming(domain.getId(), "Total Patch Smooth");
		}
		ghost_filler->fillGhost(u);
		for (const PatchInfo<D> &pinfo : domain.getPatchInfoVector()) {
			if (domain.hasTimer()) {
				domain.getTimer()->startPatchTiming(pinfo.id, domain.getId(), "Single Patch Solve");
			}
			PatchView<const double, D> f_view = f.getPatchView(pinfo.local_index);
			PatchView<double, D>       u_view = u.getPatchView(pinfo.local_index);
			solveSinglePatch(pinfo, f_view, u_view);
			if (domain.hasTimer()) {
				domain.getTimer()->stopPatchTiming(pinfo.id, domain.getId(), "Single Patch Solve");
			}
		}
		if (domain.hasTimer()) {
			domain.getTimer()->stopDomainTiming(domain.getId(), "Total Patch Smooth");
		}
	}
};
} // namespace ThunderEgg
#endif
