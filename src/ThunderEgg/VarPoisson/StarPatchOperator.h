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

#ifndef THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H
#define THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H

#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Vector.h>

namespace ThunderEgg
{
namespace VarPoisson
{
/**
 * @brief Implements a variable coefficient Laplacian f=Div[h*Grad[u]]
 *
 * h is a cell-centered coefficient
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class StarPatchOperator : public PatchOperator<D>
{
	protected:
	Vector<D> coeffs;

	constexpr int addValue(int axis) const
	{
		return (axis == 0) ? 0 : 1;
	}

	public:
	/**
	 * @brief Construct a new StarPatchOperator object
	 *
	 * @param coeffs the cell centered coefficients
	 * @param domain the Domain associated with the operator
	 * @param ghost_filler the GhostFiller to use before calling applySinglePatch
	 */
	StarPatchOperator(const Vector<D> &coeffs, const Domain<D> &domain, const GhostFiller<D> &ghost_filler)
	: PatchOperator<D>(domain, ghost_filler),
	  coeffs(coeffs)
	{
		if (domain.getNumGhostCells() < 1) {
			throw RuntimeError("StarPatchOperator needs at least one set of ghost cells");
		}
		ghost_filler.fillGhost(this->coeffs);
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
	void applySinglePatch(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view, const PatchView<double, D> &f_view) const override
	{
		PatchView<const double, D> c  = coeffs.getPatchView(pinfo.local_index);
		std::array<double, D>      h2 = pinfo.spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}
		loop<0, D - 1>([&](int axis) {
			int stride   = u_view.getStrides()[axis];
			int c_stride = c.getStrides()[axis];
			loop_over_interior_indexes<D + 1>(u_view, [&](std::array<int, D + 1> coord) {
				const double *ptr     = &u_view[coord];
				const double *c_ptr   = &c[coord];
				double        lower   = *(ptr - stride);
				double        mid     = *ptr;
				double        upper   = *(ptr + stride);
				double        c_lower = *(c_ptr - c_stride);
				double        c_mid   = *c_ptr;
				double        c_upper = *(c_ptr + c_stride);
				f_view[coord]
				= addValue(axis) * f_view[coord] + ((c_upper + c_mid) * (upper - mid) - (c_lower + c_mid) * (mid - lower)) / (2 * h2[axis]);
			});
		});
	}

	void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view) const override
	{
		for (int axis = 0; axis < D; axis++) {
			Side<D> lower_side(axis * 2);
			Side<D> upper_side(axis * 2 + 1);
			if (!pinfo.hasNbr(lower_side)) {
				View<double, D>       lower = u_view.getGhostSliceOn(lower_side, {0});
				View<const double, D> mid   = u_view.getSliceOn(lower_side, {0});
				loop_over_interior_indexes<D>(mid, [&](std::array<int, D> coord) { lower[coord] = -mid[coord]; });
			}
			if (!pinfo.hasNbr(upper_side)) {
				View<double, D>       upper = u_view.getGhostSliceOn(upper_side, {0});
				View<const double, D> mid   = u_view.getSliceOn(upper_side, {0});
				loop_over_interior_indexes<D>(mid, [&](std::array<int, D> coord) { upper[coord] = -mid[coord]; });
			}
		}
	}

	void enforceInternalBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &u_view) const override
	{
		for (int axis = 0; axis < D; axis++) {
			Side<D> lower_side(axis * 2);
			Side<D> upper_side(axis * 2 + 1);
			if (pinfo.hasNbr(lower_side)) {
				View<double, D>       lower = u_view.getGhostSliceOn(lower_side, {0});
				View<const double, D> mid   = u_view.getSliceOn(lower_side, {0});
				loop_over_interior_indexes<D>(mid, [&](std::array<int, D> coord) { lower[coord] = -mid[coord]; });
			}
			if (pinfo.hasNbr(upper_side)) {
				View<double, D>       upper = u_view.getGhostSliceOn(upper_side, {0});
				View<const double, D> mid   = u_view.getSliceOn(upper_side, {0});
				loop_over_interior_indexes<D>(mid, [&](std::array<int, D> coord) { upper[coord] = -mid[coord]; });
			}
		}
	}

	void modifyRHSForInternalBoundaryConditions(const PatchInfo<D> &              pinfo,
	                                            const PatchView<const double, D> &u_view,
	                                            const PatchView<double, D> &      f_view) const override
	{
		PatchView<const double, D> c = coeffs.getPatchView(pinfo.local_index);
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo.hasNbr(s)) {
				double                h2      = pow(pinfo.spacings[s.getAxisIndex()], 2);
				View<double, D>       fner    = f_view.getSliceOn(s, {0});
				View<double, D>       u_ghost = u_view.getGhostSliceOn(s, {0});
				View<const double, D> uner    = u_view.getSliceOn(s, {0});
				View<const double, D> c_ghost = c.getSliceOn(s, {-1});
				View<const double, D> cner    = c.getSliceOn(s, {0});
				loop_over_interior_indexes<D>(fner, [&](const std::array<int, D> &coord) {
					fner[coord] -= (u_ghost[coord] + uner[coord]) * (cner[coord] + c_ghost[coord]) / (2 * h2);
				});
			}
		}
	}
	/**
	 * @brief Helper function for adding Dirichlet boundary conditions to right hand side.
	 *
	 * @param f the right hand side vector
	 * @param gfunc the exact solution
	 * @param hfunc the coefficients
	 */
	void addDrichletBCToRHS(Vector<D> &                                          f,
	                        std::function<double(const std::array<double, D> &)> gfunc,
	                        std::function<double(const std::array<double, D> &)> hfunc)
	{
		for (int i = 0; i < f.getNumLocalPatches(); i++) {
			ComponentView<double, D> f_ld  = f.getComponentView(0, i);
			auto                     pinfo = this->getDomain().getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo.hasNbr(s)) {
					double              h2 = pow(pinfo.spacings[s.getAxisIndex()], 2);
					View<double, D - 1> ld = f_ld.getSliceOn(s, {0});
					nested_loop<D - 1>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
						std::array<double, D> real_coord;
						DomainTools::GetRealCoordBound<D>(pinfo, coord, s, real_coord);
						std::array<double, D> other_real_coord = real_coord;
						if (s.isLowerOnAxis()) {
							other_real_coord[s.getAxisIndex()] -= pinfo.spacings[s.getAxisIndex()];
						} else {
							other_real_coord[s.getAxisIndex()] += pinfo.spacings[s.getAxisIndex()];
						}
						ld[coord] -= 2 * gfunc(real_coord) * hfunc(real_coord) / h2;
					});
				}
			}
		}
	}
};
extern template class StarPatchOperator<2>;
extern template class StarPatchOperator<3>;

} // namespace VarPoisson
} // namespace ThunderEgg
#endif