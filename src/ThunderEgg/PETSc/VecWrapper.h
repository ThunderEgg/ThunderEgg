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
#include <ThunderEgg/PETSc/VecViewManager.h>
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
	 * @brief striding to next patch
	 */
	int patch_stride;
	/**
	 * @brief striding to next component
	 */
	int component_stride;
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
	std::array<int, D + 1> strides;
	/**
	 * @brief The number of (non-ghost) cells in each direction of the patch
	 */
	std::array<int, D + 1> lengths;

	/**
	 * @brief Get the Num Local Patches in this vector
	 *
	 * @param vec the Petsc vector
	 * @param lengths the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @return int the number of local patches
	 */
	static int GetNumLocalPatches(Vec vec, const std::array<int, D + 1> &lengths, int num_ghost_cells)
	{
		int patch_stride = lengths[D];
		for (size_t i = 0; i < D; i++) {
			patch_stride *= (lengths[i] + 2 * num_ghost_cells);
		}
		int vec_size;
		VecGetLocalSize(vec, &vec_size);
		return vec_size / patch_stride;
	}

	/**
	 * @brief Get the number of local cells in this vector
	 *
	 * @param vec the vector
	 * @param lengths the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @return int the number of local cells
	 */
	static int GetNumLocalCells(Vec vec, const std::array<int, D + 1> &lengths, int num_ghost_cells)
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
	VecWrapper(Vec vec, const std::array<int, D + 1> &lengths, int num_ghost_cells, bool own)
	: Vector<D>(MPI_COMM_WORLD, lengths[D], GetNumLocalPatches(vec, lengths, num_ghost_cells), GetNumLocalCells(vec, lengths, num_ghost_cells)),
	  vec(vec),
	  num_ghost_cells(num_ghost_cells),
	  own(own),
	  lengths(lengths)
	{
		strides[0] = 1;
		for (size_t i = 1; i < D; i++) {
			strides[i] = (this->lengths[i - 1] + 2 * num_ghost_cells) * strides[i - 1];
		}
		strides[D]   = strides[D - 1] * (lengths[D - 1] + 2 * num_ghost_cells);
		patch_stride = strides[D] * lengths[D];
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
	static std::shared_ptr<VecWrapper<D>> GetNewVector(std::shared_ptr<const Domain<D>> domain, int num_components)
	{
		Vec                    u;
		std::array<int, D + 1> lengths;
		for (int i = 0; i < D; i++) {
			lengths[i] = domain->getNs()[i];
		}
		lengths[D] = num_components;
		VecCreateMPI(MPI_COMM_WORLD, domain->getNumLocalCellsWithGhost() * num_components, PETSC_DETERMINE, &u);
		return std::make_shared<VecWrapper<D>>(u, lengths, domain->getNumGhostCells(), true);
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
	PatchView<double, D> getPatchView(int patch_local_index) override
	{
		std::shared_ptr<VecViewManager> ldm(new VecViewManager(vec, false));
		double *                        data = ldm->getVecView() + patch_stride * patch_local_index;
		return PatchView<double, D>(data, strides, lengths, num_ghost_cells, ldm);
	}

	PatchView<const double, D> getPatchView(int patch_local_index) const override
	{
		std::shared_ptr<VecViewManager> ldm(new VecViewManager(vec, true));
		double *                        data = ldm->getVecView() + patch_stride * patch_local_index;
		return PatchView<const double, D>(data, strides, lengths, num_ghost_cells, std::move(ldm));
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
