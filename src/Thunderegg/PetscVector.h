/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_PETSCVECTOR_H
#define THUNDEREGG_PETSCVECTOR_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/Vector.h>
#include <petscvec.h>

namespace Thunderegg
{
/**
 * @brief Manages access to petsc raw pointer data
 */
class PetscLDM : public LocalDataManager
{
	private:
	// PW<Vec> vec;
	Vec     vec;
	double *vec_view;
	bool    is_read_only;

	public:
	PetscLDM(Vec vec, bool is_read_only)
	{
		this->vec          = vec;
		this->is_read_only = is_read_only;
		if (is_read_only) {
			VecGetArrayRead(vec, const_cast<const double **>(&vec_view));
		} else {
			int state;
			VecLockGet(vec, &state);
			if (state != 0) {
				throw std::runtime_error("Vector is in state " + std::to_string(state));
			}
			VecGetArray(vec, &vec_view);
		}
	}
	~PetscLDM()
	{
		if (is_read_only) {
			VecRestoreArrayRead(vec, const_cast<const double **>(&vec_view));
		} else {
			VecRestoreArray(vec, &vec_view);
		}
	}
	double *getVecView() const
	{
		return vec_view;
	}
};
/**
 * @brief Wrap a petsc vector for use in Thunderegg
 *
 * @tparam D the number of cartesian dimensions on a patch.
 */
template <size_t D> class PetscVector : public Vector<D>

{
	private:
	int                patch_stride;
	int                num_ghost_cells;
	bool               own;
	std::array<int, D> strides;
	std::array<int, D> lengths;
	int                first_offset;
	// std::shared_ptr<void> user_data;

	public:
	PetscVector(Vec vec, int num_local_patches, const std::array<int, D> &lengths,
	            int num_ghost_cells, bool own = true)
	: Vector<D>(MPI_COMM_WORLD)
	{
		this->own             = own;
		this->lengths         = lengths;
		this->vec             = vec;
		this->num_ghost_cells = num_ghost_cells;

		strides[0]   = 1;
		first_offset = num_ghost_cells;
		for (size_t i = 1; i < D; i++) {
			strides[i] = (this->lengths[i - 1] + 2 * num_ghost_cells) * strides[i - 1];
			first_offset += strides[i] * num_ghost_cells;
		}
		patch_stride = strides[D - 1] * (lengths[D - 1] + 2 * num_ghost_cells);
		if (num_local_patches == -1) {
			VecGetLocalSize(vec, &this->num_local_patches);
			this->num_local_patches /= patch_stride;
		} else {
			this->num_local_patches = num_local_patches;
		}
		int num_cells_in_patch = 1;
		for (size_t i = 0; i < D; i++) {
			num_cells_in_patch *= this->lengths[i];
		}
		this->num_local_cells = this->num_local_patches * num_cells_in_patch;
	}

	/**
	 * @brief The petsc vector object
	 */
	Vec                 vec;
	std::vector<double> ghost_data;

	static std::shared_ptr<PetscVector<D>> GetNewVector(std::shared_ptr<const Domain<D>> domain)
	{
		Vec u;
		VecCreateMPI(MPI_COMM_WORLD, domain->getNumLocalCellsWithGhost(), PETSC_DETERMINE, &u);
		return std::shared_ptr<PetscVector<D>>(new PetscVector<D>(
		u, domain->getNumLocalPatches(), domain->getNs(), domain->getNumGhostCells()));
	}
	static std::shared_ptr<PetscVector<D - 1>> GetNewBCVector(std::shared_ptr<Domain<D>> domain)
	{
		Vec u;
		VecCreateMPI(MPI_COMM_WORLD, domain->getNumLocalBCCells(), PETSC_DETERMINE, &u);
		std::array<int, D - 1> ns;
		for (size_t i = 0; i < ns.size(); i++) {
			ns[i] = domain->getNs()[i];
		}
		return std::shared_ptr<PetscVector<D - 1>>(new PetscVector<D - 1>(u, -1, ns, 0));
	}
	/*
	static std::shared_ptr<PetscVector<D>> GetNewSchurVector(std::shared_ptr<SchurHelper<D + 1>> sh)
	{
	    Vec u;
	    VecCreateMPI(MPI_COMM_WORLD, sh->getSchurVecLocalSize(), PETSC_DETERMINE, &u);
	    return std::shared_ptr<PetscVector<D>>(new PetscVector<D>(u, sh->getLengths()));
	}
	*/

	/**
	 * @brief Destroy the Petsc Vector object
	 */
	~PetscVector()
	{
		if (own) {
			VecDestroy(&vec);
		}
	}
	/**
	 * @brief Get the LocalData object for the given local patch index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return LocalData<D> the local data object
	 */
	LocalData<D> getLocalData(int patch_local_index)
	{
		std::shared_ptr<PetscLDM> ldm(new PetscLDM(vec, false));
		double *data = ldm->getVecView() + patch_stride * patch_local_index + first_offset;
		return LocalData<D>(data, strides, lengths, num_ghost_cells, ldm);
	}
	/**
	 * @brief Get the LocalData object for the given local patch index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return LocalData<D> the local data object
	 */
	const LocalData<D> getLocalData(int patch_local_index) const
	{
		std::shared_ptr<PetscLDM> ldm(new PetscLDM(vec, true));
		double *data = ldm->getVecView() + patch_stride * patch_local_index + first_offset;
		return LocalData<D>(data, strides, lengths, num_ghost_cells, std::move(ldm));
	}
	void setNumGhostPatches(int num_ghost_patches)
	{
		ghost_data.resize(patch_stride * num_ghost_patches);
	}
	/*
	void set(double alpha)
	{
	    VecSet(vec, alpha);
	}
	void scale(double alpha)
	{
	    VecScale(vec, alpha);
	}
	void shift(double delta)
	{
	    VecShift(vec, delta);
	}
	void add(std::shared_ptr<const Vector<D>> b)
	{
	    const PetscVector<D> *b_vec = dynamic_cast<const PetscVector<D> *>(b.get());
	    if (b_vec != nullptr) {
	        VecAXPY(vec, 1, b_vec->vec);
	    } else {
	        Vector<D>::add(b);
	    }
	}
	void scaleThenAdd(double alpha, std::shared_ptr<const Vector<D>> b)
	{
	    const PetscVector<D> *b_vec = dynamic_cast<const PetscVector<D> *>(b.get());
	    if (b_vec != nullptr) {
	        VecAYPX(vec, alpha, b_vec->vec);
	    } else {
	        Vector<D>::scaleThenAdd(alpha, b);
	    }
	}
	*/
	int getNumGhostCells() const
	{
		return num_ghost_cells;
	}
};

} // namespace Thunderegg
#endif
