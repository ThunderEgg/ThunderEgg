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
 * @tparam D the number of cartesian dimensions
 */
template <size_t D> class BiCGStabPatchSolver : public PatchSolver<D>
{
	private:
	/**
	 * @brief Generates a vector that is only the size of a patch
	 */
	class SingleVG : public VectorGenerator<D>
	{
		private:
		/**
		 * @brief number of cells along each axis
		 */
		std::array<int, D> lengths;
		/**
		 * @brief number of ghost cells
		 */
		int num_ghost_cells;

		public:
		/**
		 * @brief Construct a new SingleVG object for a given patch
		 *
		 * @param pinfo the PatchInfo object for the patch
		 */
		SingleVG(std::shared_ptr<const PatchInfo<D>> pinfo)
		{
			lengths         = pinfo->ns;
			num_ghost_cells = pinfo->num_ghost_cells;
		}
		/**
		 * @brief Get a newly allocated vector
		 *
		 * @return std::shared_ptr<Vector<D>> the vector
		 */
		std::shared_ptr<Vector<D>> getNewVector()
		{
			return std::shared_ptr<Vector<D>>(new ValVector<D>(lengths, 1, num_ghost_cells));
		}
	};
	/**
	 * @brief Wrapper that provides a vector interface for the LocalData of a patch
	 */
	class SinglePatchVec : public Vector<D>
	{
		private:
		/**
		 * @brief The LocalData for the patch
		 */
		LocalData<D> ld;

		public:
		/**
		 * @brief Construct a new SinglePatchVec object
		 *
		 * @param ld the localdata for the patch
		 */
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
	/**
	 * @brief This wraps a PatchOperator object so that it only applies the operator on a specified
	 * patch
	 */
	class SinglePatchOp : public Operator<D>
	{
		private:
		/**
		 * @brief the PatchOperator that is being wrapped
		 */
		std::shared_ptr<const PatchOperator<D>> op;
		/**
		 * @brief the PatchInfo object for the patch
		 */
		std::shared_ptr<PatchInfo<D>> pinfo;

		public:
		/**
		 * @brief Construct a new SinglePatchOp object
		 *
		 * @param pinfo the patch that we want to operate on
		 * @param op the operator
		 */
		SinglePatchOp(std::shared_ptr<const PatchInfo<D>>     pinfo,
		              std::shared_ptr<const PatchOperator<D>> op)
		{
			this->pinfo = std::make_shared<PatchInfo<D>>(*pinfo);
			this->pinfo->nbr_info.fill(nullptr);
			this->op    = op;
		}
		void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
		{
			op->applySinglePatch(pinfo, x->getLocalData(0), b->getLocalData(0));
		}
	};

	/**
	 * @brief The operator for the solve
	 */
	std::shared_ptr<const PatchOperator<D>> op;

	/**
	 * @brief maximum number of iterators for each solve
	 */
	int max_it;
	/**
	 * @brief tolerance for each solve
	 */
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
        auto u_copy = vg->getNewVector();
		f_copy->copy(f_single);
        u_copy->copy(u_single);
		op->addGhostToRHS(pinfo, u, f_copy->getLocalData(0));

        // printf("Calling BiCG patch solver\n");
		BiCGStab<D>::solve(vg, single_op, u_single, f_copy, nullptr, max_it, tol);
	}
	/**
	 * @brief Generator for created a new solver for each new level in GMG
	 */
	class Generator
	{
		private:
		/**
		 * @brief A set of ghost fillers that have been generated
		 */
		std::function<std::shared_ptr<const GhostFiller<D>>(
		std::shared_ptr<const GMG::Level<D>> level)>
		filler_gen;
		/**
		 * @brief A set of PatchOperators that have been generated
		 */
		std::function<std::shared_ptr<const PatchOperator<D>>(
		std::shared_ptr<const GMG::Level<D>> level)>
		op_gen;
		/**
		 * @brief A set of PatchSolvers that have been generated
		 */
		std::map<std::shared_ptr<const Domain<D>>, std::shared_ptr<const BiCGStabPatchSolver<D>>>
		generated_solvers;
		/**
		 * @brief The tolerance
		 */
		double tol;
		/**
		 * @brief The maximum number of iterations
		 */
		int max_it;
		/**
		 * @brief timer from finest solver
		 */
		std::shared_ptr<Timer> timer;

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
			this->timer                       = solver->getTimer();
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
			if (solver != nullptr) {
				return solver;
			}
			solver.reset(new BiCGStabPatchSolver<D>(level->getDomain(), filler_gen(level),
			                                        op_gen(level), tol, max_it));
			solver->setTimer(timer);
			return solver;
		}
	};
};
} // namespace Thunderegg
#endif
