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

#include <ThunderEgg/BreakdownError.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/Vector.h>
#include <bitset>
#include <map>

namespace ThunderEgg
{
namespace Iterative
{
/**
 * @brief Solves the patches u_viewing a Iterative iterative solver on each patch
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D> class PatchSolver : public ThunderEgg::PatchSolver<D>
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
		/**
		 * @brief The number of components for each cell
		 */
		int num_components;

		public:
		/**
		 * @brief Construct a new SingleVG object for a given patch
		 *
		 * @param pinfo the PatchInfo object for the patch
		 * @param num_components the number of components for each cell
		 */
		SingleVG(const PatchInfo<D> &pinfo, int num_components)
		: lengths(pinfo.ns),
		  num_ghost_cells(pinfo.num_ghost_cells),
		  num_components(num_components)
		{
		}
		/**
		 * @brief Get a newly allocated vector
		 *
		 * @return std::shared_ptr<Vector<D>> the vector
		 */
		std::shared_ptr<Vector<D>> getNewVector() const override
		{
			return std::shared_ptr<Vector<D>>(new Vector<D>(Communicator(MPI_COMM_SELF), lengths, num_components, 1, num_ghost_cells));
		}
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
		const PatchInfo<D> &pinfo;

		public:
		/**
		 * @brief Construct a new SinglePatchOp object
		 *
		 * @param pinfo the patch that we want to operate on
		 * @param op the operator
		 */
		SinglePatchOp(const PatchInfo<D> &pinfo, std::shared_ptr<const PatchOperator<D>> op) : op(op), pinfo(pinfo) {}
		void apply(const Vector<D> &x, Vector<D> &b) const
		{
			PatchView<const double, D> x_view = x.getPatchView(0);
			PatchView<double, D>       b_view = b.getPatchView(0);
			op->enforceBoundaryConditions(pinfo, x_view);
			op->enforceZeroDirichletAtInternalBoundaries(pinfo, x_view);
			op->applySinglePatch(pinfo, x_view, b_view);
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
	 * @param op_in the PatchOperator to u_viewe
	 * @param tol_in the tolerance to u_viewe for patch solves
	 * @param max_it_in the maximum number of iterations to u_viewe for patch solves
	 * @param continue_on_breakdown continue on breakdown exception
	 */
	PatchSolver(std::shared_ptr<const Iterative::Solver<D>> solver, std::shared_ptr<const PatchOperator<D>> op_in, bool continue_on_breakdown = false)
	: ThunderEgg::PatchSolver<D>(op_in->getDomain(), op_in->getGhostFiller()),
	  solver(solver),
	  op(op_in),
	  continue_on_breakdown(continue_on_breakdown)
	{
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
		                         this->getDomain()->getNumGhostCells());

		std::array<int, D + 1> u_lengths;
		for (int i = 0; i < D + 1; i++) {
			u_lengths[i] = u_view.getEnd()[i] + 1;
		}
		Vector<D> u_single(
		Communicator(MPI_COMM_SELF), {&u_view[u_view.getGhostStart()]}, u_view.getStrides(), u_lengths, this->getDomain()->getNumGhostCells());

		Vector<D> f_copy = f_single.getZeroClone();
		f_copy.copy(f_single);
		op->modifyRHSForZeroDirichletAtInternalBoundaries(pinfo, u_view, f_copy.getPatchView(0));

		int iterations = 0;
		try {
			iterations = solver->solve(single_op, u_single, f_copy);
		} catch (const BreakdownError &err) {
			if (!continue_on_breakdown) {
				throw err;
			}
		}
		if (this->getDomain()->hasTimer()) {
			this->getDomain()->getTimer()->addIntInfo("Iterations", iterations);
		}
	}
};
} // namespace Iterative
} // namespace ThunderEgg
extern template class ThunderEgg::Iterative::PatchSolver<2>;
extern template class ThunderEgg::Iterative::PatchSolver<3>;
#endif
