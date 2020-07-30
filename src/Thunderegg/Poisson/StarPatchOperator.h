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

#ifndef THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H
#define THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H

#include <Thunderegg/DomainTools.h>
#include <Thunderegg/GMG/Level.h>
#include <Thunderegg/GhostFiller.h>
#include <Thunderegg/PatchOperator.h>
#include <Thunderegg/ValVector.h>

namespace Thunderegg
{
namespace Poisson
{
/**
 * @brief Exception that the StarPatchOperator class trows
 */
struct StarPatchOperatorException : std::runtime_error {
	StarPatchOperatorException(std::string message) : std::runtime_error(message){};
};
/**
 * @brief Implements 2nd order laplacian operator
 *
 * Supports both Dirichlet and Neumann boundary conditions
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class StarPatchOperator : public PatchOperator<D>
{
	private:
	constexpr int addValue(int axis)
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
			throw StarPatchOperatorException(
			"StarPatchOperator needs at least one set of ghost cells");
		}
	}
	void applySinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, const LocalData<D> u,
	                      LocalData<D> f) const override
	{
		std::array<double, D> h2 = pinfo->spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}

		loop<0, D - 1>([&](int axis) {
			Side<D> lower_side = Side<D>::LowerSideOnAxis(axis);
			Side<D> upper_side = Side<D>::HigherSideOnAxis(axis);
			if (!pinfo->hasNbr(lower_side)) {
				LocalData<D - 1>       lower = u.getGhostSliceOnSide(lower_side, 1);
				const LocalData<D - 1> mid   = u.getSliceOnSide(lower_side);
				if (neumann) {
					nested_loop<D - 1>(
					mid.getStart(), mid.getEnd(),
					[&](std::array<int, D - 1> coord) { lower[coord] = mid[coord]; });
				} else {
					nested_loop<D - 1>(
					mid.getStart(), mid.getEnd(),
					[&](std::array<int, D - 1> coord) { lower[coord] = -mid[coord]; });
				}
			}
			if (!pinfo->hasNbr(upper_side)) {
				LocalData<D - 1>       upper = u.getGhostSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid   = u.getSliceOnSide(upper_side);
				if (neumann) {
					nested_loop<D - 1>(
					mid.getStart(), mid.getEnd(),
					[&](std::array<int, D - 1> coord) { upper[coord] = mid[coord]; });
				} else {
					nested_loop<D - 1>(
					mid.getStart(), mid.getEnd(),
					[&](std::array<int, D - 1> coord) { upper[coord] = -mid[coord]; });
				}
			}
			int stride = u.getStrides()[axis];
			nested_loop<D>(u.getStart(), u.getEnd(), [&](std::array<int, D> coord) {
				const double *ptr   = u.getPtr(coord);
				double        lower = *(ptr - stride);
				double        mid   = *ptr;
				double        upper = *(ptr + stride);
				f[coord] = addValue(axis) * f[coord] + (upper - 2 * mid + lower) / h2[axis];
			});
		});
	}
	void addGhostToRHS(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                   LocalData<D> f) const
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				double           h2      = pow(pinfo->spacings[s.getAxisIndex()], 2);
				LocalData<D - 1> f_inner = f.getSliceOnSide(s);
				LocalData<D - 1> u_ghost = u.getSliceOnSide(s, -1);
				LocalData<D - 1> u_inner = u.getSliceOnSide(s);
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
			LocalData<D> f_ld  = f->getLocalData(i);
			auto         pinfo = this->domain->getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo->hasNbr(s)) {
					double           h2 = pow(pinfo->spacings[s.getAxisIndex()], 2);
					LocalData<D - 1> ld = f_ld.getSliceOnSide(s);
					nested_loop<D - 1>(
					ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
						std::array<double, D> real_coord;
						DomainTools<D>::getRealCoordBound(pinfo, coord, s, real_coord);
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
			LocalData<D> f_ld  = f->getLocalData(i);
			auto         pinfo = this->domain->getPatchInfoVector()[i];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo->hasNbr(s)) {
					double           h  = pinfo->spacings[s.getAxisIndex()];
					LocalData<D - 1> ld = f_ld.getSliceOnSide(s);
					if (s.isLowerOnAxis()) {
						nested_loop<D - 1>(
						ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
							std::array<double, D> real_coord;
							DomainTools<D>::getRealCoordBound(pinfo, coord, s, real_coord);
							ld[coord] += gfunc_grad[s.getAxisIndex()](real_coord) / h;
						});
					} else {
						nested_loop<D - 1>(
						ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
							std::array<double, D> real_coord;
							DomainTools<D>::getRealCoordBound(pinfo, coord, s, real_coord);
							ld[coord] -= gfunc_grad[s.getAxisIndex()](real_coord) / h;
						});
					}
				}
			}
		}
	}
	/**
	 * @brief Generator for GMG levels.
	 *
	 * Will use same interpolation scheme for coefficions as the
	 * interpolator in GMG.
	 */
	class Generator
	{
		private:
		/**
		 * @brief generator for ghost fillers
		 */
		std::function<std::shared_ptr<const GhostFiller<D>>(
		std::shared_ptr<const GMG::Level<D>> level)>
		filler_gen;
		/**
		 * @brief Generated operators are stored here.
		 */
		std::map<std::shared_ptr<const Domain<D>>, std::shared_ptr<const StarPatchOperator<D>>>
		generated_operators;

		public:
		/**
		 * @brief Construct a new StarPatchOperator generator
		 *
		 * @param finest_op the finest star pach operator
		 * @param filler_gen returns a GhostFiller for a given level
		 */
		Generator(std::shared_ptr<const StarPatchOperator<D>> finest_op,
		          std::function<
		          std::shared_ptr<const GhostFiller<D>>(std::shared_ptr<const GMG::Level<D>> level)>
		          filler_gen_in)
		: filler_gen(filler_gen_in)
		{
			generated_operators[finest_op->domain] = finest_op;
		}
		/**
		 * @brief Return a StarPatchOperator for a given level
		 *
		 * @param level the level in GMG
		 * @return std::shared_ptr<const StarPatchOperator<D>> the operator
		 */
		std::shared_ptr<const StarPatchOperator<D>>
		operator()(std::shared_ptr<const GMG::Level<D>> level)
		{
			auto &coarser_op = generated_operators[level->getDomain()];
			if (coarser_op != nullptr) {
				return coarser_op;
			}

			std::shared_ptr<const Domain<D>> finer_domain = level->getFiner()->getDomain();
			auto                             finer_op     = generated_operators[finer_domain];
			coarser_op.reset(new StarPatchOperator<D>(level->getDomain(), filler_gen(level)));
			return coarser_op;
		}
	};
};
extern template class StarPatchOperator<2>;
extern template class StarPatchOperator<3>;

} // namespace Poisson
} // namespace Thunderegg
#endif