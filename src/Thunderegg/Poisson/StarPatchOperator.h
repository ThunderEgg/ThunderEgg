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
template <size_t D> class StarPatchOperator : public PatchOperator<D>
{
	private:
	constexpr int addValue(int axis)
	{
		return (axis == 0) ? 0 : 1;
	}

	public:
	StarPatchOperator(std::shared_ptr<const Domain<D>>      domain_in,
	                  std::shared_ptr<const GhostFiller<D>> ghost_filler_in)
	: PatchOperator<D>(domain_in, ghost_filler_in)
	{
		if (this->domain->getNumGhostCells() < 1) {
			throw 88;
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
			Side<D> lower_side(axis * 2);
			Side<D> upper_side(axis * 2 + 1);
			if (!pinfo->hasNbr(lower_side)) {
				LocalData<D - 1>       lower = u.getGhostSliceOnSide(lower_side, 1);
				const LocalData<D - 1> mid   = u.getSliceOnSide(lower_side);
				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					lower[coord] = -mid[coord];
				});
			}
			if (!pinfo->hasNbr(upper_side)) {
				LocalData<D - 1>       upper = u.getGhostSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid   = u.getSliceOnSide(upper_side);
				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					upper[coord] = -mid[coord];
				});
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
	 * @param domain The domain associated with the vector
	 * @param f the right hand side vector
	 * @param gfunc the exact solution
	 * @param hfunc the coefficients
	 */
	static void addDrichletBCToRHS(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> f,
	                               std::function<double(const std::array<double, D> &)> gfunc)
	{
		for (int i = 0; i < f->getNumLocalPatches(); i++) {
			LocalData<D> f_ld  = f->getLocalData(i);
			auto         pinfo = domain->getPatchInfoVector()[i];
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
		          filler_gen)
		{
			generated_operators[finest_op->domain] = finest_op;
			this->filler_gen                       = filler_gen;
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