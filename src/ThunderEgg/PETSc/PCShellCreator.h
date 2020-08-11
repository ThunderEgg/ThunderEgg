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

#ifndef THUNDEREGG_PETSC_PCSHELLCREATOR_H
#define THUNDEREGG_PETSC_PCSHELLCREATOR_H
#include <ThunderEgg/PETSc/MatShellCreator.h>
#include <petscpc.h>
namespace ThunderEgg
{
namespace PETSc
{
/**
 * @brief Wraps a ThunderEgg operator for use as a PETSc PC
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class PCShellCreator
{
	private:
	/**
	 * @brief The operator we are wrapping
	 */
	std::shared_ptr<Operator<D>> op;
	/**
	 * @brief The associated VectorGenerator
	 */
	std::shared_ptr<VectorGenerator<D>> vg;
	/**
	 * @brief The Mat associated with the preconditioner operator
	 */
	Mat A;

	/**
	 * @brief Construct a new PCShellCreator object
	 *
	 * @param op_in the Operator we are wrapping
	 * @param vg_in the VectorGenerator we are wrapping
	 * @param A_in the Mat associated with the preconditioner Operator
	 */
	PCShellCreator(std::shared_ptr<Operator<D>> op_in, std::shared_ptr<VectorGenerator<D>> vg_in,
	               Mat A_in)
	: op(op_in), vg(vg_in), A(A_in)
	{
	}
	PCShellCreator(const PCShellCreator &) = delete;
	PCShellCreator &operator=(const PCShellCreator &) = delete;
	PCShellCreator(PCShellCreator &&) noexcept        = delete;
	PCShellCreator &operator=(PCShellCreator &&) noexcept = delete;
	/**
	 * @brief Destroy the PCShellCreator object
	 */
	~PCShellCreator()
	{
		MatDestroy(&A);
	}
	/**
	 * @brief Apply the preconditioner
	 *
	 * @param A the preconditioner
	 * @param x the input vector
	 * @param b the output vector
	 * @return int PETSc error code
	 */
	static int applyPC(PC A, Vec x, Vec b)
	{
		PCShellCreator<D> *psc = nullptr;
		PCShellGetContext(A, (void **) &psc);

		auto te_x = psc->vg->getNewVector();
		auto te_b = psc->vg->getNewVector();

		// petsc vectors don't have gost padding for patchs, so this is neccesary
		const double *x_view;
		VecGetArrayRead(x, &x_view);
		int index = 0;
		for (int p_index = 0; p_index < te_x->getNumLocalPatches(); p_index++) {
			LocalData<D> ld = te_x->getLocalData(p_index);
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
				ld[coord] = x_view[index];
				index++;
			});
		}

		VecRestoreArrayRead(x, &x_view);

		psc->op->apply(te_x, te_b);

		double *b_view;
		VecGetArray(b, &b_view);
		index = 0;
		for (int p_index = 0; p_index < te_b->getNumLocalPatches(); p_index++) {
			const LocalData<D> ld = te_b->getLocalData(p_index);
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
				b_view[index] = ld[coord];
				index++;
			});
		}
		VecRestoreArray(b, &b_view);

		return 0;
	}
	/**
	 * @brief deallocate the wrapper
	 *
	 * @param A the wrapped object
	 * @return int PETSc error
	 */
	static int destroyPC(PC A)
	{
		PCShellCreator<D> *psc = nullptr;
		PCShellGetContext(A, (void **) &psc);
		delete psc;
		return 0;
	}

	public:
	/**
	 * @brief Construct a new PCShell object
	 *
	 * @param prec the preconditioner that is being wrapped
	 * @param op the operator that is being preconditioned
	 * @param vg the VectorGenerator associated with the operators
	 * @return PC the wrapped PC, you are responsible for calling PCDestroy on this object
	 */
	static PC GetNewPCShell(std::shared_ptr<Operator<D>> prec, std::shared_ptr<Operator<D>> op,
	                        std::shared_ptr<VectorGenerator<D>> vg)
	{
		Mat                A   = MatShellCreator<D>::GetNewMatShell(op, vg);
		PCShellCreator<D> *psc = new PCShellCreator(prec, vg, A);
		PC                 P;
		PCCreate(MPI_COMM_WORLD, &P);
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, psc);
		PCShellSetApply(P, applyPC);
		PCShellSetDestroy(P, destroyPC);
		PCSetOperators(P, A, A);
		return P;
	}
};
} // namespace PETSc
} // namespace ThunderEgg
#endif
