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

#ifndef THUNDEREGG_SCHUR_BICGSTABSOLVER_H
#define THUNDEREGG_SCHUR_BICGSTABSOLVER_H

#include <ThunderEgg/BiCGStab.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/ValVector.h>
#include <bitset>
#include <map>

namespace ThunderEgg
{
/**
 * @brief Solves the patches using a BiCGStab iterative solver on each patch
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D> class BiCGStabPatchSolver : public PatchSolver<D>
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
		std::shared_ptr<Vector<D>> getNewVector() const override
		{
			return std::shared_ptr<Vector<D>>(
			new ValVector<D>(MPI_COMM_SELF, lengths, 1, 1, num_ghost_cells));
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
		SinglePatchVec(const LocalData<D> &ld)
		: Vector<D>(MPI_COMM_SELF, 1, 1, GetNumLocalCells(ld)), ld(ld)
		{
		}
		LocalData<D> getLocalData(int component_index, int local_patch_id) override
		{
			return ld;
		}
		const LocalData<D> getLocalData(int component_index, int local_patch_id) const override
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
			this->op = op;
		}
		void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
		{
			auto xs = x->getLocalDatas(0);
			auto bs = b->getLocalDatas(0);
			op->applySinglePatch(pinfo, xs, bs);
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
	 * @brief Construct a new BiCGStabPatchSolver object
	 *
	 * @param op_in the PatchOperator to use
	 * @param tol_in the tolerance to use for patch solves
	 * @param max_it_in the maximum number of iterations to use for patch solves
	 */
	BiCGStabPatchSolver(std::shared_ptr<const PatchOperator<D>> op_in, double tol_in = 1e-12,
	                    int max_it_in = 1000)
	: PatchSolver<D>(op_in->getDomain(), op_in->getGhostFiller()), op(op_in), max_it(max_it_in),
	  tol(tol_in)
	{
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
		const std::vector<LocalData<D>> us = {u};
		std::vector<LocalData<D>>       fs = {f_copy->getLocalData(0, 0)};
		op->addGhostToRHS(pinfo, us, fs);
		// u_single->set(0);

		// printf("Calling BiCG patch solver\n");
		BiCGStab<D>::solve(vg, single_op, u_single, f_copy, nullptr, max_it, tol);
	}
};
} // namespace ThunderEgg
#endif
