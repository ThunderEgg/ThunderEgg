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

#ifndef THUNDEREGG_PETSC_MATSHELLCREATOR_H
#define THUNDEREGG_PETSC_MATSHELLCREATOR_H
#include <ThunderEgg/Operator.h>
#include <ThunderEgg/VectorGenerator.h>
#include <petscmat.h>
namespace ThunderEgg
{
namespace PETSc
{
template <size_t D> class MatShellCreator
{
	private:
	std::shared_ptr<Operator<D>>        op;
	std::shared_ptr<VectorGenerator<D>> vg;

	MatShellCreator(std::shared_ptr<Operator<D>> op_in, std::shared_ptr<VectorGenerator<D>> vg_in)
	: op(op_in), vg(vg_in)
	{
	}
	static int applyMat(Mat A, Vec x, Vec b)
	{
		MatShellCreator<D> *msc = nullptr;
		MatShellGetContext(A, &msc);

		auto te_x = msc->vg->getNewVector();
		auto te_b = msc->vg->getNewVector();

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

		msc->op->apply(te_x, te_b);

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
	static int destroyMat(Mat A)
	{
		MatShellCreator<D> *msc = nullptr;
		MatShellGetContext(A, &msc);
		delete msc;
		return 0;
	}

	public:
	static Mat GetNewMatShell(std::shared_ptr<Operator<D>>        op,
	                          std::shared_ptr<VectorGenerator<D>> vg)
	{
		MatShellCreator<D> *msc = new MatShellCreator(op, vg);
		int                 m   = vg->getNewVector()->getNumLocalCells();
		Mat                 A;
		MatCreateShell(MPI_COMM_WORLD, m, m, PETSC_DETERMINE, PETSC_DETERMINE, msc, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) applyMat);
		MatShellSetOperation(A, MATOP_DESTROY, (void (*)(void)) destroyMat);
		return A;
	}
};
} // namespace PETSc
} // namespace ThunderEgg
#endif
