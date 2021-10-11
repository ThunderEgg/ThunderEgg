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

#ifndef THUNDEREGG_POISSON_STARPATCHOPERATOR_H
#define THUNDEREGG_POISSON_STARPATCHOPERATOR_H
/**
 * @file
 *
 * @brief StarPatchOperator class
 */

#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Vector.h>

namespace ThunderEgg::Poisson
{
/**
 * @brief Implements 2nd order finite-difference Laplacian operator
 *
 * Supports both Dirichlet and Neumann boundary conditions
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class StarPatchOperator : public PatchOperator<D>
{
	private:
	constexpr int addValue(int axis) const
	{
		return (axis == 0) ? 0 : 1;
	}
	bool neumann;

	public:
	/**
	 * @brief Construct a new StarPatchOperator object
	 *
	 * @param domain the Domain that the operator is associated with
	 * @param ghost_filler the GhostFiller to use before calling applySinglePatch
	 * @param neumann whether or not to use Neumann boundary conditions
	 */
	StarPatchOperator(const Domain<D> &domain, const GhostFiller<D> &ghost_filler, bool neumann = false)
	: PatchOperator<D>(domain, ghost_filler),
	  neumann(neumann)
	{
		if (this->getDomain().getNumGhostCells() < 1) {
			throw RuntimeError("StarPatchOperator needs at least one set of ghost cells");
		}
	}
	/**
	 * @brief Get a clone of this operator
	 *
	 * @return StarPatchOperator<D>* a newly allocated copy of this operator
	 */
	StarPatchOperator<D> *clone() const override
	{
		return new StarPatchOperator<D>(*this);
	}
	void
	applySinglePatch(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view, const PatchView<double, D> &f_view, bool internal) const
	{
		if (internal) {
			enforceInternalBoundaryConditions(pinfo, u_view);
		}

		enforceBoundaryConditions(pinfo, u_view);

		std::array<double, D> h2 = pinfo.spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}

		Loop::Unroll<0, D - 1>([&](int axis) {
			int stride = u_view.getStrides()[axis];
			Loop::OverInteriorIndexes<D + 1>(u_view, [&](std::array<int, D + 1> coord) {
				const double *ptr   = &u_view[coord];
				double        lower = *(ptr - stride);
				double        mid   = *ptr;
				double        upper = *(ptr + stride);
				f_view[coord]       = addValue(axis) * f_view[coord] + (upper - 2 * mid + lower) / h2[axis];
			});
		});
	}
	void applySinglePatch(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view, const PatchView<double, D> &f_view) const override
	{
		applySinglePatch(pinfo, u_view, f_view, false);
	}
	void applySinglePatchWithInternalBoundaryConditions(const PatchInfo<D>               &pinfo,
	                                                    const PatchView<const double, D> &u_view,
	                                                    const PatchView<double, D>       &f_view) const override
	{
		applySinglePatch(pinfo, u_view, f_view, true);
	}
	void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view) const
	{
		for (int axis = 0; axis < D; axis++) {
			Side<D> lower_side = LowerSideOnAxis<D>(axis);
			Side<D> upper_side = HigherSideOnAxis<D>(axis);
			if (!pinfo.hasNbr(lower_side)) {
				View<double, D>       lower     = u_view.getGhostSliceOn(lower_side, {0});
				View<const double, D> lower_mid = u_view.getSliceOn(lower_side, {0});
				if (neumann) {
					Loop::OverInteriorIndexes<D>(lower_mid, [&](std::array<int, D> coord) { lower[coord] = lower_mid[coord]; });
				} else {
					Loop::OverInteriorIndexes<D>(lower_mid, [&](std::array<int, D> coord) { lower[coord] = -lower_mid[coord]; });
				}
			}
			if (!pinfo.hasNbr(upper_side)) {
				View<double, D>       upper     = u_view.getGhostSliceOn(upper_side, {0});
				View<const double, D> upper_mid = u_view.getSliceOn(upper_side, {0});
				if (neumann) {
					Loop::OverInteriorIndexes<D>(upper_mid, [&](std::array<int, D> coord) { upper[coord] = upper_mid[coord]; });
				} else {
					Loop::OverInteriorIndexes<D>(upper_mid, [&](std::array<int, D> coord) { upper[coord] = -upper_mid[coord]; });
				}
			}
		}
	}
	void enforceInternalBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view) const
	{
		for (int axis = 0; axis < D; axis++) {
			Side<D> lower_side = LowerSideOnAxis<D>(axis);
			Side<D> upper_side = HigherSideOnAxis<D>(axis);
			if (pinfo.hasNbr(lower_side)) {
				View<double, D>       lower     = u_view.getGhostSliceOn(lower_side, {0});
				View<const double, D> lower_mid = u_view.getSliceOn(lower_side, {0});
				Loop::OverInteriorIndexes<D>(lower_mid, [&](std::array<int, D> coord) { lower[coord] = -lower_mid[coord]; });
			}
			if (pinfo.hasNbr(upper_side)) {
				View<double, D>       upper     = u_view.getGhostSliceOn(upper_side, {0});
				View<const double, D> upper_mid = u_view.getSliceOn(upper_side, {0});
				Loop::OverInteriorIndexes<D>(upper_mid, [&](std::array<int, D> coord) { upper[coord] = -upper_mid[coord]; });
			}
		}
	}
	void modifyRHSForInternalBoundaryConditions(const PatchInfo<D>               &pinfo,
	                                            const PatchView<const double, D> &u_view,
	                                            const PatchView<double, D>       &f_view) const override
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo.hasNbr(s)) {
				double                h2      = pow(pinfo.spacings[s.getAxisIndex()], 2);
				View<double, D>       f_inner = f_view.getSliceOn(s, {0});
				View<const double, D> u_ghost = u_view.getSliceOn(s, {-1});
				View<const double, D> u_inner = u_view.getSliceOn(s, {0});
				Loop::OverInteriorIndexes<D>(f_inner,
				                             [&](const std::array<int, D> &coord) { f_inner[coord] -= (u_ghost[coord] + u_inner[coord]) / h2; });
			}
		}
	}
	/**
	 * @brief Helper function for adding Dirichlet boundary conditions to right hand side.
	 *
	 * @param f the right hand side vector
	 * @param gfunc the exact solution
	 */
	void addDrichletBCToRHS(Vector<D> &f, std::function<double(const std::array<double, D> &)> gfunc)
	{
		for (int i = 0; i < f.getNumLocalPatches(); i++) {
			ComponentView<double, D> f_ld  = f.getComponentView(0, i);
			auto                     pinfo = this->getDomain().getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo.hasNbr(s)) {
					double              h2 = pow(pinfo.spacings[s.getAxisIndex()], 2);
					View<double, D - 1> ld = f_ld.getSliceOn(s, {0});
					Loop::Nested<D - 1>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
						std::array<double, D> real_coord;
						DomainTools::GetRealCoordBound<D>(pinfo, coord, s, real_coord);
						ld[coord] += -2.0 * gfunc(real_coord) / h2;
					});
				}
			}
		}
	}
	/**
	 * @brief Helper function for adding Neumann boundary conditions to right hand side.
	 *
	 * @param f the right hand side vector
	 * @param gfunc the exact solution
	 * @param gfunc_grad the gradient of gfunc
	 */
	void addNeumannBCToRHS(Vector<D>                                                          &f,
	                       std::function<double(const std::array<double, D> &)>                gfunc,
	                       std::array<std::function<double(const std::array<double, D> &)>, D> gfunc_grad)
	{
		for (int i = 0; i < f.getNumLocalPatches(); i++) {
			ComponentView<double, D> f_ld  = f.getComponentView(0, i);
			auto                     pinfo = this->getDomain().getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo.hasNbr(s)) {
					double              h  = pinfo.spacings[s.getAxisIndex()];
					View<double, D - 1> ld = f_ld.getSliceOn(s, {0});
					if (s.isLowerOnAxis()) {
						Loop::Nested<D - 1>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
							std::array<double, D> real_coord;
							DomainTools::GetRealCoordBound<D>(pinfo, coord, s, real_coord);
							ld[coord] += gfunc_grad[s.getAxisIndex()](real_coord) / h;
						});
					} else {
						Loop::Nested<D - 1>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
							std::array<double, D> real_coord;
							DomainTools::GetRealCoordBound<D>(pinfo, coord, s, real_coord);
							ld[coord] -= gfunc_grad[s.getAxisIndex()](real_coord) / h;
						});
					}
				}
			}
		}
	}
};
extern template class StarPatchOperator<2>;
extern template class StarPatchOperator<3>;
} // namespace ThunderEgg::Poisson
#endif