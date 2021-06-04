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

#ifndef THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H
#define THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H

#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/ValVector.h>

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
	std::shared_ptr<const Vector<D>> coeffs;

	constexpr int addValue(int axis) const
	{
		return (axis == 0) ? 0 : 1;
	}

	public:
	/**
	 * @brief Construct a new StarPatchOperator object
	 *
	 * @param coeffs_in the cell centered coefficients
	 * @param domain_in the Domain associated with the operator
	 * @param ghost_filler_in the GhostFiller to use before calling applySinglePatch
	 */
	StarPatchOperator(std::shared_ptr<const Vector<D>>      coeffs_in,
	                  std::shared_ptr<const Domain<D>>      domain_in,
	                  std::shared_ptr<const GhostFiller<D>> ghost_filler_in)
	: PatchOperator<D>(domain_in, ghost_filler_in),
	  coeffs(coeffs_in)
	{
		if (this->domain->getNumGhostCells() < 1) {
			throw RuntimeError("StarPatchOperator needs at least one set of ghost cells");
		}
		this->ghost_filler->fillGhost(this->coeffs);
	}
	void applySinglePatch(const PatchInfo<D> &pinfo, const std::vector<ComponentView<D>> &us, std::vector<ComponentView<D>> &fs) const override
	{
		const ComponentView<D> c  = coeffs->getComponentView(0, pinfo.local_index);
		std::array<double, D>  h2 = pinfo.spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}
		loop<0, D - 1>([&](int axis) {
			int stride   = us[0].getStrides()[axis];
			int c_stride = c.getStrides()[axis];
			nested_loop<D>(us[0].getStart(), us[0].getEnd(), [&](std::array<int, D> coord) {
				const double *ptr     = &us[0][coord];
				const double *c_ptr   = &c[coord];
				double        lower   = *(ptr - stride);
				double        mid     = *ptr;
				double        upper   = *(ptr + stride);
				double        c_lower = *(c_ptr - c_stride);
				double        c_mid   = *c_ptr;
				double        c_upper = *(c_ptr + c_stride);
				fs[0][coord]
				= addValue(axis) * fs[0][coord] + ((c_upper + c_mid) * (upper - mid) - (c_lower + c_mid) * (mid - lower)) / (2 * h2[axis]);
			});
		});
	}

	void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const std::vector<ComponentView<D>> &us) const override
	{
		for (int axis = 0; axis < D; axis++) {
			Side<D> lower_side(axis * 2);
			Side<D> upper_side(axis * 2 + 1);
			if (!pinfo.hasNbr(lower_side)) {
				ComponentView<D - 1>       lower = us[0].getSliceOn(lower_side, {-1});
				const ComponentView<D - 1> mid   = us[0].getSliceOn(lower_side, {0});
				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) { lower[coord] = -mid[coord]; });
			}
			if (!pinfo.hasNbr(upper_side)) {
				ComponentView<D - 1>       upper = us[0].getSliceOn(upper_side, {-1});
				const ComponentView<D - 1> mid   = us[0].getSliceOn(upper_side, {0});
				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) { upper[coord] = -mid[coord]; });
			}
		}
	}

	void enforceZeroDirichletAtInternalBoundaries(const PatchInfo<D> &pinfo, const std::vector<ComponentView<D>> &us) const override
	{
		for (int axis = 0; axis < D; axis++) {
			Side<D> lower_side(axis * 2);
			Side<D> upper_side(axis * 2 + 1);
			if (pinfo.hasNbr(lower_side)) {
				ComponentView<D - 1>       lower = us[0].getSliceOn(lower_side, {-1});
				const ComponentView<D - 1> mid   = us[0].getSliceOn(lower_side, {0});
				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) { lower[coord] = -mid[coord]; });
			}
			if (pinfo.hasNbr(upper_side)) {
				ComponentView<D - 1>       upper = us[0].getSliceOn(upper_side, {-1});
				const ComponentView<D - 1> mid   = us[0].getSliceOn(upper_side, {0});
				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) { upper[coord] = -mid[coord]; });
			}
		}
	}

	void modifyRHSForZeroDirichletAtInternalBoundaries(const PatchInfo<D> &                 pinfo,
	                                                   const std::vector<ComponentView<D>> &us,
	                                                   std::vector<ComponentView<D>> &      fs) const override
	{
		const ComponentView<D> c = coeffs->getComponentView(0, pinfo.local_index);
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo.hasNbr(s)) {
				double                     h2      = pow(pinfo.spacings[s.getAxisIndex()], 2);
				ComponentView<D - 1>       f_inner = fs[0].getSliceOn(s, {0});
				ComponentView<D - 1>       u_ghost = us[0].getSliceOn(s, {-1});
				const ComponentView<D - 1> u_inner = us[0].getSliceOn(s, {0});
				const ComponentView<D - 1> c_ghost = c.getSliceOn(s, {-1});
				const ComponentView<D - 1> c_inner = c.getSliceOn(s, {0});
				nested_loop<D - 1>(f_inner.getStart(), f_inner.getEnd(), [&](const std::array<int, D - 1> &coord) {
					f_inner[coord] -= (u_ghost[coord] + u_inner[coord]) * (c_inner[coord] + c_ghost[coord]) / (2 * h2);
					u_ghost[coord] = 0;
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
	void addDrichletBCToRHS(std::shared_ptr<Vector<D>>                           f,
	                        std::function<double(const std::array<double, D> &)> gfunc,
	                        std::function<double(const std::array<double, D> &)> hfunc)
	{
		for (int i = 0; i < f->getNumLocalPatches(); i++) {
			ComponentView<D> f_ld  = f->getComponentView(0, i);
			auto             pinfo = this->domain->getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo.hasNbr(s)) {
					double               h2 = pow(pinfo.spacings[s.getAxisIndex()], 2);
					ComponentView<D - 1> ld = f_ld.getSliceOn(s, {0});
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