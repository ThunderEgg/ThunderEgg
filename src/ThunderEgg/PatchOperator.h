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
 *  This program is distributed in the hope that it will be u_vieweful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef THUNDEREGG_PATCHOPERATOR_H
#define THUNDEREGG_PATCHOPERATOR_H
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/Operator.h>
#include <ThunderEgg/Vector.h>
namespace ThunderEgg
{
/**
 * @brief This is an Operator where derived classes only have to implement the two virtual functions
 * that operate on single patch.
 *
 * @tparam D the number of Cartesian dimensions.
 */
template <int D> class PatchOperator : public Operator<D>
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
	 * @brief Construct a new Patch Operator object
	 *
	 *  This sets the Domain and GhostFiller
	 *
	 * @param domain  the Domain
	 * @param ghost_filler the GhostFiller
	 */
	PatchOperator(const Domain<D> &domain, const GhostFiller<D> &ghost_filler) : domain(domain), ghost_filler(ghost_filler.clone()) {}
	/**
	 * @brief Clone this patch operator
	 *
	 * @return PatchOperator<D>* a newly allocated copy of this patch operator
	 */
	virtual PatchOperator<D> *clone() const override = 0;
	/**
	 * @brief Destroy the PatchOperator object
	 */
	virtual ~PatchOperator() {}

	/**
	 * @brief Apply the operator to a single patch
	 *
	 * The ghost values in u will be updated to the latest values, and should not need to be modified
	 *
	 * @param pinfo  the patch
	 * @param u_view the right hand side
	 * @param f_view the left hand side
	 * @param treat_interior_boundary_as_dirichlet if true, the stencil of the patch should be
	 * modified so that the interior boundaries are assumed to be zero, and the ghost values should
	 * not be u_viewed
	 */
	virtual void applySinglePatch(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view, const PatchView<double, D> &f_view) const = 0;

	/**
	 * @brief modify values in ghost cells in order to enforce boundary conditions
	 *
	 * @param pinfo the patch info
	 * @param u_view the left hand side
	 */
	virtual void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view) const = 0;

	/**
	 * @brief modify values in ghost cells order to enforce the internal boundary conditions
	 *
	 * @param pinfo the patch info
	 * @param u_view the left hand side
	 */
	virtual void enforceInternalBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view) const = 0;

	/**
	 * @brief Treat the internal patch boundaries as domain boundaires and modify
	 * RHS accordingly.
	 *
	 * This will be u_viewed in patch solvers to formulate a RHS for the individual patch to solve for.
	 *
	 * @param pinfo the patch
	 * @param u_view the left hand side
	 * @param f_view the right hand side
	 */
	virtual void modifyRHSForInternalBoundaryConditions(const PatchInfo<D> &              pinfo,
	                                                    const PatchView<const double, D> &u_view,
	                                                    const PatchView<double, D> &      f_view) const = 0;

	/**
	 * @brief Apply the operator
	 *
	 * This will update the ghost values in u, and then will call applySinglePatch for each patch
	 *
	 * @param u the left hand side
	 * @param f the right hand side
	 */
	void apply(const Vector<D> &u, Vector<D> &f) const override
	{
		if constexpr (ENABLE_DEBUG) {
			if (u.getNumLocalPatches() != this->domain->getNumLocalPatches()) {
				throw RuntimeError("u vector is incorrect length");
			}
			if (f.getNumLocalPatches() != this->domain->getNumLocalPatches()) {
				throw RuntimeError("f vector is incorrect length");
			}
		}
		f.setWithGhost(0);
		ghost_filler->fillGhost(u);
		for (const PatchInfo<D> &pinfo : domain.getPatchInfoVector()) {
			PatchView<const double, D> u_view = u.getPatchView(pinfo.local_index);
			enforceBoundaryConditions(pinfo, u_view);
			PatchView<double, D> f_view = f.getPatchView(pinfo.local_index);
			applySinglePatch(pinfo, u_view, f_view);
		}
	}
	/**
	 * @brief Get the Domain object associated with this PatchOperator
	 */
	const Domain<D> &getDomain() const
	{
		return domain;
	}
	/**
	 * @brief Get the GhostFiller object associated with this PatchOperator
	 */
	const GhostFiller<D> &getGhostFiller() const
	{
		return *ghost_filler;
	}
};
} // namespace ThunderEgg
extern template class ThunderEgg::PatchOperator<2>;
extern template class ThunderEgg::PatchOperator<3>;
#endif
