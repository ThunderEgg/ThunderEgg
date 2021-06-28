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
#include <petscmat.h>
namespace ThunderEgg
{
namespace PETSc
{
/**
 * @brief Wraps a ThunderEgg operator for use as a PETSc Mat
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class MatShellCreator
{
	private:
	/**
	 * @brief the operator we are wrapping
	 */
	std::shared_ptr<const Operator<D>> op;

	/**
	 * @brief get a new vector
	 */
	std::function<Vector<D>()> getNewVector;

	/**
	 * @brief Construct a new MatShellCreator object
	 *
	 * @param op the Operator
	 */
	explicit MatShellCreator(std::shared_ptr<const Operator<D>> op, const std::function<Vector<D>()> &vector_allocator)
	: op(op),
	  getNewVector(vector_allocator)
	{
	}
	/**
	 * @brief Apply the PETSC MatShell
	 *
	 * @param A the MatShell
	 * @param x the x vector (input)
	 * @param b the b vector (output)
	 * @return int PETSC error
	 */
	static int applyMat(Mat A, Vec x, Vec b)
	{
		MatShellCreator<D> *msc = nullptr;
		MatShellGetContext(A, &msc);

		Vector<2> te_x = msc->getNewVector();
		Vector<2> te_b = msc->getNewVector();

		// petsc vectors don't have gost padding for patchs, so this is neccesary
		const double *x_view;
		VecGetArrayRead(x, &x_view);
		int index = 0;
		for (int p_index = 0; p_index < te_x.getNumLocalPatches(); p_index++) {
			for (int c = 0; c < te_x.getNumComponents(); c++) {
				ComponentView<double, D> ld = te_x.getComponentView(c, p_index);
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
					ld[coord] = x_view[index];
					index++;
				});
			}
		}

		VecRestoreArrayRead(x, &x_view);

		msc->op->apply(te_x, te_b);

		double *b_view;
		VecGetArray(b, &b_view);
		index = 0;
		for (int p_index = 0; p_index < te_b.getNumLocalPatches(); p_index++) {
			for (int c = 0; c < te_b.getNumComponents(); c++) {
				const ComponentView<const double, D> ld = te_b.getComponentView(c, p_index);
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
					b_view[index] = ld[coord];
					index++;
				});
			}
		}
		VecRestoreArray(b, &b_view);

		return 0;
	}
	/**
	 * @brief Deallocate the wrapper
	 *
	 * @param A the wrapped operator
	 * @return int PETSc error code
	 */
	static int destroyMat(Mat A)
	{
		MatShellCreator<D> *msc = nullptr;
		MatShellGetContext(A, &msc);
		delete msc;
		return 0;
	}

	public:
	/**
	 * @brief Get a new MatShell for use with PETSc
	 *
	 * @param op the operator we are wrapping
	 * @param vg the associated vectorGenerator for that object
	 * @return Mat the wrapped operator, you are responsible for calling MatDestroy on this
	 */
	static Mat GetNewMatShell(std::shared_ptr<const Operator<D>> op, const std::function<Vector<D>()> &vector_allocator)
	{
		MatShellCreator<D> *msc = new MatShellCreator(op, vector_allocator);
		Vector<D>           vec = vector_allocator();
		int                 m   = vec.getNumLocalCells() * vec.getNumComponents();
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
