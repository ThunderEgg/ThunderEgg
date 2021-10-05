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

#ifndef THUNDEREGG_ITERATIVE_PATCHSOLVER_H
#define THUNDEREGG_ITERATIVE_PATCHSOLVER_H
/**
 * @file
 *
 * @brief PatchSolver class
 */

#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/Iterative/BreakdownError.h>
#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/Vector.h>
#include <bitset>
#include <map>

namespace ThunderEgg::Iterative
{
/**
 * @brief Solves the patches using an iterative Solver on each patch
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D> class PatchSolver : public ThunderEgg::PatchSolver<D>
{
	private:
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
		const PatchInfo<D> &pinfo;

		public:
		/**
		 * @brief Construct a new SinglePatchOp object
		 *
		 * @param pinfo the patch that we want to operate on
		 * @param op the operator
		 */
		SinglePatchOp(const PatchInfo<D> &pinfo, std::shared_ptr<const PatchOperator<D>> op) : op(op), pinfo(pinfo) {}
		void apply(const Vector<D> &x, Vector<D> &b) const override
		{
			PatchView<const double, D> x_view = x.getPatchView(0);
			PatchView<double, D>       b_view = b.getPatchView(0);
			op->applySinglePatchWithInternalBoundaryConditions(pinfo, x_view, b_view);
		}
		/**
		 * @brief Get a clone of this operator
		 *
		 * @return SinglePatchOp* a newly allocated copy
		 */
		SinglePatchOp *clone() const override
		{
			return new SinglePatchOp(*this);
		}
	};

	/**
	 * @brief The iterative solver being u_viewed
	 */
	std::shared_ptr<const Solver<D>> solver;

	/**
	 * @brief The operator for the solve
	 */
	std::shared_ptr<const PatchOperator<D>> op;

	/**
	 * @brief whether or not to continue on BreakDownError
	 */
	bool continue_on_breakdown;

	public:
	/**
	 * @brief Construct a new IterativePatchSolver object
	 *
	 * @param tol_in the tolerance to use for patch solves
	 * @param max_it_in the maximum number of iterations to use for patch solves
	 */

	/**
	 * @brief Construct a new Patch Solver object
	 *
	 * @param solver the solver to use
	 * @param op the PatchOperator to use
	 * @param continue_on_breakdown continue on breakdown exception
	 */
	PatchSolver(const Iterative::Solver<D> &solver, const PatchOperator<D> &op, bool continue_on_breakdown = false)
	: ThunderEgg::PatchSolver<D>(op.getDomain(), op.getGhostFiller()),
	  solver(solver.clone()),
	  op(op.clone()),
	  continue_on_breakdown(continue_on_breakdown)
	{
	}
	/**
	 * @brief Clone this patch solver
	 *
	 * @return PatchSolver<D>* a newly allocated copy of this patch solver
	 */
	PatchSolver<D> *clone() const override
	{
		return new PatchSolver<D>(*this);
	}
	void solveSinglePatch(const PatchInfo<D> &pinfo, const PatchView<const double, D> &f_view, const PatchView<double, D> &u_view) const override
	{
		SinglePatchOp single_op(pinfo, op);

		std::array<int, D + 1> f_lengths;
		for (int i = 0; i < D + 1; i++) {
			f_lengths[i] = f_view.getEnd()[i] + 1;
		}
		const Vector<D> f_single(Communicator(MPI_COMM_SELF),
		                         {const_cast<double *>(&f_view[f_view.getGhostStart()])},
		                         f_view.getStrides(),
		                         f_lengths,
		                         this->getDomain().getNumGhostCells());

		std::array<int, D + 1> u_lengths;
		for (int i = 0; i < D + 1; i++) {
			u_lengths[i] = u_view.getEnd()[i] + 1;
		}
		Vector<D> u_single(
		Communicator(MPI_COMM_SELF), {&u_view[u_view.getGhostStart()]}, u_view.getStrides(), u_lengths, this->getDomain().getNumGhostCells());

		Vector<D> f_copy = f_single.getZeroClone();
		f_copy.copy(f_single);
		op->modifyRHSForInternalBoundaryConditions(pinfo, u_view, f_copy.getPatchView(0));

		int iterations = 0;
		try {
			iterations = solver->solve(single_op, u_single, f_copy);
		} catch (const BreakdownError &err) {
			if (!continue_on_breakdown) {
				throw err;
			}
		}
		if (this->getDomain().hasTimer()) {
			this->getDomain().getTimer()->addIntInfo("Iterations", iterations);
		}
	}
};
} // namespace ThunderEgg::Iterative
extern template class ThunderEgg::Iterative::PatchSolver<2>;
extern template class ThunderEgg::Iterative::PatchSolver<3>;
#endif
