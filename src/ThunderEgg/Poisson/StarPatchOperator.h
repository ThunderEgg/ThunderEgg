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

#ifndef THUNDEREGG_POISSON_STARPATCHOPERATOR_H
#define THUNDEREGG_POISSON_STARPATCHOPERATOR_H

#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/ValVector.h>

namespace ThunderEgg
{
namespace Poisson
{
/**
 * @brief Implements 2nd order laplacian operator
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
	 * @param domain_in the Domain that the operator is associated with
	 * @param ghost_filler_in the GhostFiller to use before calling applySinglePatch
	 * @param neumann_in whether or not to use Neumann boundary conditions
	 */
	StarPatchOperator(std::shared_ptr<const Domain<D>>      domain_in,
	                  std::shared_ptr<const GhostFiller<D>> ghost_filler_in,
	                  bool                                  neumann_in = false)
	: PatchOperator<D>(domain_in, ghost_filler_in), neumann(neumann_in)
	{
		if (this->domain->getNumGhostCells() < 1) {
			throw RuntimeError("StarPatchOperator needs at least one set of ghost cells");
		}
	}
	void applySinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                      const std::vector<LocalData<D>> &us, std::vector<LocalData<D>> &fs,
	                      bool treat_interior_boundary_as_dirichlet) const override
	{
		std::array<double, D> h2 = pinfo->spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}

		loop<0, D - 1>([&](int axis) {
			Side<D>                lower_side = Side<D>::LowerSideOnAxis(axis);
			Side<D>                upper_side = Side<D>::HigherSideOnAxis(axis);
			LocalData<D - 1>       lower      = us[0].getGhostSliceOnSide(lower_side, 1);
			const LocalData<D - 1> lower_mid  = us[0].getSliceOnSide(lower_side);
			if (!pinfo->hasNbr(lower_side) && neumann) {
				nested_loop<D - 1>(
				lower_mid.getStart(), lower_mid.getEnd(),
				[&](std::array<int, D - 1> coord) { lower[coord] = lower_mid[coord]; });
			} else if (!pinfo->hasNbr(lower_side) || treat_interior_boundary_as_dirichlet) {
				nested_loop<D - 1>(
				lower_mid.getStart(), lower_mid.getEnd(),
				[&](std::array<int, D - 1> coord) { lower[coord] = -lower_mid[coord]; });
			}
			LocalData<D - 1>       upper     = us[0].getGhostSliceOnSide(upper_side, 1);
			const LocalData<D - 1> upper_mid = us[0].getSliceOnSide(upper_side);
			if (!pinfo->hasNbr(upper_side) && neumann) {
				nested_loop<D - 1>(
				upper_mid.getStart(), upper_mid.getEnd(),
				[&](std::array<int, D - 1> coord) { upper[coord] = upper_mid[coord]; });
			} else if (!pinfo->hasNbr(upper_side) || treat_interior_boundary_as_dirichlet) {
				nested_loop<D - 1>(
				upper_mid.getStart(), upper_mid.getEnd(),
				[&](std::array<int, D - 1> coord) { upper[coord] = -upper_mid[coord]; });
			}
			int stride = us[0].getStrides()[axis];
			nested_loop<D>(us[0].getStart(), us[0].getEnd(), [&](std::array<int, D> coord) {
				const double *ptr   = us[0].getPtr(coord);
				double        lower = *(ptr - stride);
				double        mid   = *ptr;
				double        upper = *(ptr + stride);
				fs[0][coord] = addValue(axis) * fs[0][coord] + (upper - 2 * mid + lower) / h2[axis];
			});
		});
	}
	void addGhostToRHS(std::shared_ptr<const PatchInfo<D>> pinfo,
	                   const std::vector<LocalData<D>> &   us,
	                   std::vector<LocalData<D>> &         fs) const override
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				double                 h2      = pow(pinfo->spacings[s.getAxisIndex()], 2);
				LocalData<D - 1>       f_inner = fs[0].getSliceOnSide(s);
				LocalData<D - 1>       u_ghost = us[0].getSliceOnSide(s, -1);
				const LocalData<D - 1> u_inner = us[0].getSliceOnSide(s);
				nested_loop<D - 1>(f_inner.getStart(), f_inner.getEnd(),
				                   [&](const std::array<int, D - 1> &coord) {
					                   f_inner[coord] -= (u_ghost[coord] + u_inner[coord]) / h2;
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
	 */
	void addDrichletBCToRHS(std::shared_ptr<Vector<D>>                           f,
	                        std::function<double(const std::array<double, D> &)> gfunc)
	{
		for (int i = 0; i < f->getNumLocalPatches(); i++) {
			LocalData<D> f_ld  = f->getLocalData(0, i);
			auto         pinfo = this->domain->getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo->hasNbr(s)) {
					double           h2 = pow(pinfo->spacings[s.getAxisIndex()], 2);
					LocalData<D - 1> ld = f_ld.getSliceOnSide(s);
					nested_loop<D - 1>(
					ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
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
	void addNeumannBCToRHS(
	std::shared_ptr<Vector<D>> f, std::function<double(const std::array<double, D> &)> gfunc,
	std::array<std::function<double(const std::array<double, D> &)>, D> gfunc_grad)
	{
		for (int i = 0; i < f->getNumLocalPatches(); i++) {
			LocalData<D> f_ld  = f->getLocalData(0, i);
			auto         pinfo = this->domain->getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo->hasNbr(s)) {
					double           h  = pinfo->spacings[s.getAxisIndex()];
					LocalData<D - 1> ld = f_ld.getSliceOnSide(s);
					if (s.isLowerOnAxis()) {
						nested_loop<D - 1>(
						ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
							std::array<double, D> real_coord;
							DomainTools::GetRealCoordBound<D>(pinfo, coord, s, real_coord);
							ld[coord] += gfunc_grad[s.getAxisIndex()](real_coord) / h;
						});
					} else {
						nested_loop<D - 1>(
						ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
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

} // namespace Poisson
} // namespace ThunderEgg
#endif