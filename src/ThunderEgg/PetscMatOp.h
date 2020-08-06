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

#ifndef THUNDEREGG_PETSCMATOP_H
#define THUNDEREGG_PETSCMATOP_H

#include <ThunderEgg/Operator.h>
#include <ThunderEgg/PetscVector.h>
#include <petscmat.h>

namespace ThunderEgg
{
/**
 * @brief Base class for operators
 */
template <size_t D> class PetscMatOp : public Operator<D>
{
	private:
	Mat A;
	Vec getPetscVecWithoutGhost(std::shared_ptr<const Vector<D>> vec) const
	{
		Vec                                   petsc_vec;
		std::shared_ptr<const PetscVector<D>> petsc_vec_ptr
		= std::dynamic_pointer_cast<const PetscVector<D>>(vec);
		if (petsc_vec_ptr != nullptr && petsc_vec_ptr->getNumGhostCells() == 0) {
			petsc_vec = petsc_vec_ptr->vec;
		} else {
			// have to create a new petsc vector without ghostcells for petsc call
			VecCreateMPI(vec->getMPIComm(), vec->getNumLocalCells(), PETSC_DETERMINE, &petsc_vec);
		}
		return petsc_vec;
	}
	void copyToPetscVec(std::shared_ptr<const Vector<D>> vec, Vec petsc_vec) const
	{
		std::shared_ptr<const PetscVector<D>> petsc_vec_ptr
		= std::dynamic_pointer_cast<const PetscVector<D>>(vec);
		if (petsc_vec_ptr == nullptr || petsc_vec_ptr->getNumGhostCells() > 0) {
			double *petsc_vec_view;
			size_t  curr_index = 0;
			VecGetArray(petsc_vec, &petsc_vec_view);
			for (int i = 0; i < vec->getNumLocalPatches(); i++) {
				const LocalData<D> ld = vec->getLocalData(i);
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
					petsc_vec_view[curr_index] = ld[coord];
					curr_index++;
				});
			}
			VecRestoreArray(petsc_vec, &petsc_vec_view);
		}
	}
	void copyToVec(Vec petsc_vec, std::shared_ptr<Vector<D>> vec) const
	{
		std::shared_ptr<PetscVector<D>> petsc_vec_ptr
		= std::dynamic_pointer_cast<PetscVector<D>>(vec);
		if (petsc_vec_ptr == nullptr || petsc_vec_ptr->getNumGhostCells() > 0) {
			const double *petsc_vec_view;
			size_t        curr_index = 0;
			VecGetArrayRead(petsc_vec, &petsc_vec_view);
			for (int i = 0; i < vec->getNumLocalPatches(); i++) {
				LocalData<D> ld = vec->getLocalData(i);
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
					ld[coord] = petsc_vec_view[curr_index];
					curr_index++;
				});
			}
			VecRestoreArrayRead(petsc_vec, &petsc_vec_view);
		}
	}
	void destroyPetscVec(std::shared_ptr<const Vector<D>> vec, Vec petsc_vec) const
	{
		std::shared_ptr<const PetscVector<D>> petsc_vec_ptr
		= std::dynamic_pointer_cast<const PetscVector<D>>(vec);
		if (petsc_vec_ptr == nullptr || petsc_vec_ptr->getNumGhostCells() > 0) {
			VecDestroy(&petsc_vec);
		}
	}

	public:
	PetscMatOp(Mat A)
	{
		this->A = A;
	}
	/**
	 * @brief Apply Petsc matrix
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
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
} // namespace ThunderEgg
#endif