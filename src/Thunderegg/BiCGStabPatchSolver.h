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

#ifndef THUNDEREGG_SCHUR_BICGSTABSOLVER_H
#define THUNDEREGG_SCHUR_BICGSTABSOLVER_H

#include <Thunderegg/BiCGStab.h>
#include <Thunderegg/Domain.h>
#include <Thunderegg/GMG/Level.h>
#include <Thunderegg/PatchOperator.h>
#include <Thunderegg/PatchSolver.h>
#include <Thunderegg/ValVector.h>
#include <bitset>
#include <map>

namespace Thunderegg
{
/**
 * @brief Solves the patches using a BiCGStab iterative solver on each patch
 *
 * @tparam D
 */
template <size_t D> class BiCGStabPatchSolver : public PatchSolver<D>
{
	private:
	class SingleVG : public VectorGenerator<D>
	{
		private:
		std::array<int, D> lengths;
		int                num_ghost_cells;

		public:
		SingleVG(std::shared_ptr<const PatchInfo<D>> pinfo)
		{
			lengths         = pinfo->ns;
			num_ghost_cells = pinfo->num_ghost_cells;
		}
		std::shared_ptr<Vector<D>> getNewVector()
		{
			return std::shared_ptr<Vector<D>>(new ValVector<D>(lengths, 1, num_ghost_cells));
		}
	};
	class SinglePatchVec : public Vector<D>
	{
		private:
		LocalData<D> ld;

		public:
		SinglePatchVec(const LocalData<D> &ld)
		{
			this->num_local_patches = 1;
			this->ld                = ld;
		}
		LocalData<D> getLocalData(int local_patch_id)
		{
			return ld;
		}
		const LocalData<D> getLocalData(int local_patch_id) const
		{
			return ld;
		}
		void setNumGhostPatches(int) {}
	};
	class SinglePatchOp : public Operator<D>
	{
		private:
		std::shared_ptr<const PatchOperator<D>> op;
		std::shared_ptr<const PatchInfo<D>>     pinfo;

		public:
		SinglePatchOp(std::shared_ptr<const PatchInfo<D>>     pinfo,
		              std::shared_ptr<const PatchOperator<D>> op)
		{
			this->pinfo = pinfo;
			this->op    = op;
		}
		void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
		{
			op->applySinglePatch(pinfo, x->getLocalData(0), b->getLocalData(0));
		}
	};
	std::shared_ptr<const PatchOperator<D>> op;

	protected:
	int    max_it;
	double tol;

	public:
	/**
	 * @brief Construct a new BiCGStabSolver object
	 *
	 * @param op the operator to use
	 * @param tol the tolerance
	 * @param max_it the maximum iterations
	 */
	BiCGStabPatchSolver(std::shared_ptr<const Domain<D>>        domain,
	                    std::shared_ptr<const GhostFiller<D>>   ghost_filler,
	                    std::shared_ptr<const PatchOperator<D>> op, double tol = 1e-12,
	                    int max_it = 1000)
	{
		this->domain       = domain;
		this->ghost_filler = ghost_filler;
		this->op           = op;
		this->tol          = tol;
		this->max_it       = max_it;
	}
	void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                      const LocalData<D> f) const override
	{
		std::shared_ptr<SinglePatchOp>      single_op(new SinglePatchOp(pinfo, op));
		std::shared_ptr<VectorGenerator<D>> vg(new SingleVG(pinfo));

		std::shared_ptr<Vector<D>> f_single(new SinglePatchVec(f));
		std::shared_ptr<Vector<D>> u_single(new SinglePatchVec(u));

		auto f_copy = vg->getNewVector();
		f_copy->copy(f_single);
		op->addGhostToRHS(pinfo, u, f_copy->getLocalData(0));

		BiCGStab<D>::solve(vg, single_op, u_single, f_copy, nullptr, max_it, tol);
	}
	class Generator
	{
		private:
		std::function<std::shared_ptr<const GhostFiller<D>>(
		std::shared_ptr<const GMG::Level<D>> level)>
		filler_gen;
		std::function<std::shared_ptr<const PatchOperator<D>>(
		std::shared_ptr<const GMG::Level<D>> level)>
		op_gen;
		std::map<std::shared_ptr<const Domain<D>>, std::shared_ptr<const BiCGStabPatchSolver<D>>>
		       generated_solvers;
		double tol;
		int    max_it;

		public:
		/**
		 * @brief Construct a new Generator object
		 *
		 * @param solver the solver for the finest level
		 * @param filler_gen GhostFiller generator
		 * @param op_gen Operator generator
		 */
		Generator(std::shared_ptr<const BiCGStabPatchSolver<D>> solver,
		          std::function<
		          std::shared_ptr<const GhostFiller<D>>(std::shared_ptr<const GMG::Level<D>> level)>
		          filler_gen,
		          std::function<std::shared_ptr<const PatchOperator<D>>(
		          std::shared_ptr<const GMG::Level<D>> level)>
		          op_gen)
		{
			this->tol                         = solver->tol;
			this->max_it                      = solver->max_it;
			this->filler_gen                  = filler_gen;
			this->op_gen                      = op_gen;
			generated_solvers[solver->domain] = solver;
		}
		/**
		 * @brief get a new BiCGStabPatchSolver for given gmg level
		 *
		 * @param level the level
		 * @return std::shared_ptr<const BiCGStabPatchSolver<D>>
		 */
		std::shared_ptr<const BiCGStabPatchSolver<D>>
		operator()(std::shared_ptr<const GMG::Level<D>> level)
		{
			auto solver = generated_solvers[level->getDomain()];
			if (solver != nullptr) { return solver; }
			solver.reset(new BiCGStabPatchSolver<D>(level->getDomain(), filler_gen(level),
			                                        op_gen(level), tol, max_it));
			return solver;
		}
	};
};
} // namespace Thunderegg
#endif
