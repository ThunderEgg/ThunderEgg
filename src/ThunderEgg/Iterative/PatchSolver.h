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

#ifndef THUNDEREGG_ITERATIVE_PATCHSOLVER_H
#define THUNDEREGG_ITERATIVE_PATCHSOLVER_H

#include <ThunderEgg/BreakdownError.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/ValVector.h>
#include <bitset>
#include <map>

namespace ThunderEgg
{
namespace Iterative
{
/**
 * @brief Solves the patches using a Iterative iterative solver on each patch
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
		SingleVG(std::shared_ptr<const PatchInfo<D>> pinfo, int num_components)
		: lengths(pinfo->ns), num_ghost_cells(pinfo->num_ghost_cells),
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
			return std::shared_ptr<Vector<D>>(
			new ValVector<D>(MPI_COMM_SELF, lengths, num_ghost_cells, num_components, 1));
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
		std::vector<LocalData<D>> lds;

		/**
		 * @brief Get the number of local cells in the LocalData object
		 *
		 * @param ld the LocalData object
		 * @return int the number of local cells
		 */
		static int GetNumLocalCells(const LocalData<D> &ld)
		{
			int patch_stride = 1;
			for (size_t i = 0; i < D; i++) {
				patch_stride *= ld.getLengths()[i];
			}
			return patch_stride;
		}

		public:
		/**
		 * @brief Construct a new SinglePatchVec object
		 *
		 * @param ld_in the localdata for the patch
		 */
		SinglePatchVec(const std::vector<LocalData<D>> &lds)
		: Vector<D>(MPI_COMM_SELF, lds.size(), 1, GetNumLocalCells(lds[0])), lds(lds)
		{
		}
		LocalData<D> getLocalData(int component_index, int local_patch_id) override
		{
			return lds[component_index];
		}
		const LocalData<D> getLocalData(int component_index, int local_patch_id) const override
		{
			return lds[component_index];
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
		std::shared_ptr<const PatchInfo<D>> pinfo;

		public:
		/**
		 * @brief Construct a new SinglePatchOp object
		 *
		 * @param pinfo the patch that we want to operate on
		 * @param op the operator
		 */
		SinglePatchOp(std::shared_ptr<const PatchInfo<D>>     pinfo,
		              std::shared_ptr<const PatchOperator<D>> op)
		: op(op), pinfo(pinfo)
		{
		}
		void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
		{
			auto xs = x->getLocalDatas(0);
			auto bs = b->getLocalDatas(0);
			op->applySinglePatch(pinfo, xs, bs, true);
		}
	};

	/**
	 * @brief The iterative solver being used
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
	 * @param op_in the PatchOperator to use
	 * @param tol_in the tolerance to use for patch solves
	 * @param max_it_in the maximum number of iterations to use for patch solves
	 * @param continue_on_breakdown continue on breakdown exception
	 */
	PatchSolver(std::shared_ptr<const Iterative::Solver<D>> solver,
	            std::shared_ptr<const PatchOperator<D>> op_in, bool continue_on_breakdown = false)
	: ThunderEgg::PatchSolver<D>(op_in->getDomain(), op_in->getGhostFiller()), solver(solver),
	  op(op_in), continue_on_breakdown(continue_on_breakdown)
	{
	}
	void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                      const std::vector<LocalData<D>> &   fs,
	                      std::vector<LocalData<D>> &         us) const override
	{
		std::shared_ptr<SinglePatchOp>      single_op(new SinglePatchOp(pinfo, op));
		std::shared_ptr<VectorGenerator<D>> vg(new SingleVG(pinfo, fs.size()));

		std::shared_ptr<Vector<D>> f_single(new SinglePatchVec(fs));
		std::shared_ptr<Vector<D>> u_single(new SinglePatchVec(us));

		auto f_copy = vg->getNewVector();
		f_copy->copy(f_single);
		auto f_copy_lds = f_copy->getLocalDatas(0);
		op->addGhostToRHS(pinfo, us, f_copy_lds);

		int iterations = 0;
		try {
			iterations = solver->solve(vg, single_op, u_single, f_copy);
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
