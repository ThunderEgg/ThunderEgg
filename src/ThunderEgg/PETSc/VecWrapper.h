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

#ifndef THUNDEREGG_PETSC_VECWRAPPER_H
#define THUNDEREGG_PETSC_VECWRAPPER_H
#include <ThunderEgg/PETSc/VecLocalDataManager.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>

namespace ThunderEgg
{
namespace PETSc
{
/**
 * @brief Wrap a PETSc Vec for use as a ThunderEgg Vector
 *
 * @tparam D the number of Cartesian dimensions on a patch.
 */
template <int D> class VecWrapper : public Vector<D>

{
	private:
	/**
	 * @brief The petsc vector object
	 */
	Vec vec;
	/**
	 * @brief The number of cells (including ghost) in each patch
	 */
	int patch_stride;
	/**
	 * @brief The number of ghost cell rows on each side of the patch
	 */
	int num_ghost_cells;
	/**
	 * @brief Whether or not to deallocate PETSc Vec on destruction
	 */
	bool own;
	/**
	 * @brief The striding for each patch
	 */
	std::array<int, D> strides;
	/**
	 * @brief The number of (non-ghost) cells in each direction of the patch
	 */
	std::array<int, D> lengths;
	/**
	 * @brief The offset to the first non-ghost cell value
	 */
	int first_offset;

	/**
	 * @brief Get the Num Local Patches in this vector
	 *
	 * @param vec the Petsc vector
	 * @param lengths the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @return int the number of local patches
	 */
	static int GetNumLocalPatches(Vec vec, const std::array<int, D> &lengths, int num_ghost_cells)
	{
		int patch_stride = 1;
		for (size_t i = 0; i < D; i++) {
			patch_stride *= (lengths[i] + 2 * num_ghost_cells);
		}
		int num_cells;
		VecGetLocalSize(vec, &num_cells);
		return num_cells / patch_stride;
	}

	/**
	 * @brief Get the number of local cells in this vector
	 *
	 * @param vec the vector
	 * @param lengths the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @return int the number of local cells
	 */
	static int GetNumLocalCells(Vec vec, const std::array<int, D> &lengths, int num_ghost_cells)
	{
		int num_cells_in_patch = 1;
		for (size_t i = 0; i < D; i++) {
			num_cells_in_patch *= lengths[i];
		}
		return num_cells_in_patch * GetNumLocalPatches(vec, lengths, num_ghost_cells);
	}

	public:
	/**
	 * @brief Construct a new VecWrapper object
	 *
	 * @param vec the allocated PETSc vector
	 * @param lengths the number of (non-ghost) cells in each direction of the patch
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @param num_components the number of components for each cell
	 * @param own whether or not to deallocate the PETSc vector when this object is destroyed.
	 */
	VecWrapper(Vec vec, const std::array<int, D> &lengths, int num_components, int num_ghost_cells,
	           bool own)
	: Vector<D>(MPI_COMM_WORLD, 1, GetNumLocalPatches(vec, lengths, num_ghost_cells),
	            GetNumLocalCells(vec, lengths, num_ghost_cells)),
	  vec(vec), num_ghost_cells(num_ghost_cells), own(own), lengths(lengths)
	{
		strides[0]   = 1;
		first_offset = num_ghost_cells;
		for (size_t i = 1; i < D; i++) {
			strides[i] = (this->lengths[i - 1] + 2 * num_ghost_cells) * strides[i - 1];
			first_offset += strides[i] * num_ghost_cells;
		}
		patch_stride = strides[D - 1] * (lengths[D - 1] + 2 * num_ghost_cells);
	}
	VecWrapper(const VecWrapper &) = delete;
	VecWrapper &operator=(const VecWrapper &) = delete;
	VecWrapper(VecWrapper &&) noexcept        = delete;
	VecWrapper &operator=(VecWrapper &&) noexcept = delete;

	/**
	 * @brief Get a new VecWrapper object for a given Domain
	 *
	 * @param domain the Domain
	 * @return std::shared_ptr<VecWrapper<D>> the resulting VecWrapper
	 */
	static std::shared_ptr<VecWrapper<D>> GetNewVector(std::shared_ptr<const Domain<D>> domain,
	                                                   int num_components)
	{
		Vec u;
		VecCreateMPI(MPI_COMM_WORLD, domain->getNumLocalCellsWithGhost(), PETSC_DETERMINE, &u);
		return std::make_shared<VecWrapper<D>>(u, domain->getNs(), num_components,
		                                       domain->getNumGhostCells(), true);
	}
	/**
	 * @brief Get a new boundary condition vector for a given Domain
	 *
	 * @param domain the Domain
	 * @return std::shared_ptr<VecWrapper<D - 1>> the resulting vector
	 */
	static std::shared_ptr<VecWrapper<D - 1>> GetNewBCVector(std::shared_ptr<Domain<D>> domain)
	{
		Vec u;
		VecCreateMPI(MPI_COMM_WORLD, domain->getNumLocalBCCells(), PETSC_DETERMINE, &u);
		std::array<int, D - 1> ns;
		for (size_t i = 0; i < ns.size(); i++) {
			ns[i] = domain->getNs()[i];
		}
		return std::make_shared<VecWrapper<D - 1>>(u, ns, 1, 0, true);
	}
	/**
	 * @brief Destroy the VecWrapper object
	 */
	~VecWrapper()
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
	LocalData<D> getLocalData(int component_index, int patch_local_index) override
	{
		std::shared_ptr<VecLocalDataManager> ldm(new VecLocalDataManager(vec, false));
		double *data = ldm->getVecView() + patch_stride * patch_local_index + first_offset;
		return LocalData<D>(data, strides, lengths, num_ghost_cells, ldm);
	}

	/**
	 * @brief Get the LocalData object for the given local patch index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return LocalData<D> the local data object
	 */
	const LocalData<D> getLocalData(int component_index, int patch_local_index) const override
	{
		std::shared_ptr<VecLocalDataManager> ldm(new VecLocalDataManager(vec, true));
		double *data = ldm->getVecView() + patch_stride * patch_local_index + first_offset;
		return LocalData<D>(data, strides, lengths, num_ghost_cells, std::move(ldm));
	}
	int getNumGhostCells() const
	{
		return num_ghost_cells;
	}
	/**
	 * @brief Get the underlying Vec object
	 *
	 * @return Vec the Vec object
	 */
	Vec getVec() const
	{
		return vec;
	}
};

} // namespace PETSc
} // namespace ThunderEgg
#endif
