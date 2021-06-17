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

#ifndef THUNDEREGG_PETSC_MATWRAPPER_H
#define THUNDEREGG_PETSC_MATWRAPPER_H

#include <ThunderEgg/Operator.h>
#include <ThunderEgg/PETSc/VecWrapper.h>
#include <petscmat.h>

namespace ThunderEgg
{
namespace PETSc
{
/**
 * @brief Wraps a PETSc Mat object for use as a ThunderEgg Operator
 */
template <int D> class MatWrapper : public Operator<D>
{
	private:
	/**
	 * @brief The PETSc matrix
	 */
	Mat A;
	/**
	 * @brief Allocate a PETSc Vec without the additional padding for ghost cells
	 *
	 * @param vec the ThunderEgg Vector
	 * @return Vec the PETSc Vec
	 */
	Vec getPetscVecWithoutGhost(const Vector<D> &vec) const
	{
		Vec                         petsc_vec;
		const PETSc::VecWrapper<D> *petsc_vec_ptr = dynamic_cast<const PETSc::VecWrapper<D> *>(&vec);
		if (petsc_vec_ptr != nullptr && petsc_vec_ptr->getNumGhostCells() == 0) {
			petsc_vec = petsc_vec_ptr->getVec();
		} else {
			// have to create a new petsc vector without ghostcells for petsc call
			VecCreateMPI(vec.getMPIComm(), vec.getNumLocalCells() * vec.getNumComponents(), PETSC_DETERMINE, &petsc_vec);
		}
		return petsc_vec;
	}
	/**
	 * @brief Copy to a PETSc Vec from a ThunderEgg Vector
	 *
	 * @param vec the ThunderEgg Vector
	 * @param petsc_vec the PETSc Vec to copy to
	 */
	void copyToPetscVec(const Vector<D> &vec, Vec petsc_vec) const
	{
		const PETSc::VecWrapper<D> *petsc_vec_ptr = dynamic_cast<const PETSc::VecWrapper<D> *>(&vec);
		if (petsc_vec_ptr == nullptr || petsc_vec_ptr->getNumGhostCells() > 0) {
			double *petsc_vec_view;
			size_t  curr_index = 0;
			VecGetArray(petsc_vec, &petsc_vec_view);
			for (int i = 0; i < vec.getNumLocalPatches(); i++) {
				for (int c = 0; c < vec.getNumComponents(); c++) {
					const ComponentView<const double, D> ld = vec.getComponentView(c, i);
					nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
						petsc_vec_view[curr_index] = ld[coord];
						curr_index++;
					});
				}
			}
			VecRestoreArray(petsc_vec, &petsc_vec_view);
		}
	}
	/**
	 * @brief Copy from a PETSc Vec to a ThunderEgg Vector
	 *
	 * @param petsc_vec the PETSc Vec
	 * @param vec the ThunderEgg Vector
	 */
	void copyToVec(Vec petsc_vec, Vector<D> &vec) const
	{
		PETSc::VecWrapper<D> *petsc_vec_ptr = dynamic_cast<PETSc::VecWrapper<D> *>(&vec);
		if (petsc_vec_ptr == nullptr || petsc_vec_ptr->getNumGhostCells() > 0) {
			const double *petsc_vec_view;
			size_t        curr_index = 0;
			VecGetArrayRead(petsc_vec, &petsc_vec_view);
			for (int i = 0; i < vec.getNumLocalPatches(); i++) {
				for (int c = 0; c < vec.getNumComponents(); c++) {
					ComponentView<double, D> ld = vec.getComponentView(c, i);
					nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
						ld[coord] = petsc_vec_view[curr_index];
						curr_index++;
					});
				}
			}
			VecRestoreArrayRead(petsc_vec, &petsc_vec_view);
		}
	}
	/**
	 * @brief Deallocate the PETSc vector
	 *
	 * @param vec the corresponding ThunderEgg Vector
	 * @param petsc_vec the PETSc Vec to deallocate
	 */
	void destroyPetscVec(const Vector<D> &vec, Vec petsc_vec) const
	{
		const PETSc::VecWrapper<D> *petsc_vec_ptr = dynamic_cast<const PETSc::VecWrapper<D> *>(&vec);
		if (petsc_vec_ptr == nullptr || petsc_vec_ptr->getNumGhostCells() > 0) {
			VecDestroy(&petsc_vec);
		}
	}

	public:
	/**
	 * @brief Construct a new MatWrapper object
	 *
	 * This object will not deallocate the PETSc Mat, you are responsible for deallocating it.
	 *
	 * @param A_in the PETSc Mat to wrap
	 */
	explicit MatWrapper(Mat A_in) : A(A_in) {}
	void apply(const Vector<D> &x, Vector<D> &b) const override
	{
		Vec petsc_x = getPetscVecWithoutGhost(x);
		copyToPetscVec(x, petsc_x);

		Vec petsc_b = getPetscVecWithoutGhost(b);

		MatMult(A, petsc_x, petsc_b);

		copyToVec(petsc_b, b);

		destroyPetscVec(x, petsc_x);
		destroyPetscVec(b, petsc_b);
	}
};
} // namespace PETSc
} // namespace ThunderEgg
#endif