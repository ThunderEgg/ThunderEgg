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
template <size_t D> class PCShellCreator
{
	private:
	std::shared_ptr<Operator<D>>        op;
	std::shared_ptr<VectorGenerator<D>> vg;
	Mat                                 A;

	PCShellCreator(std::shared_ptr<Operator<D>> op_in, std::shared_ptr<VectorGenerator<D>> vg_in,
	               Mat A_in)
	: op(op_in), vg(vg_in), A(A_in)
	{
	}
	PCShellCreator(const PCShellCreator &) = delete;
	~PCShellCreator()
	{
		MatDestroy(&A);
	}
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
	static int destroyPC(PC A)
	{
		PCShellCreator<D> *psc = nullptr;
		PCShellGetContext(A, (void **) &psc);
		delete psc;
		return 0;
	}

	public:
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
