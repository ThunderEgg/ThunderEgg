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

#ifndef THUNDEREGG_VALVECTOR_H
#define THUNDEREGG_VALVECTOR_H
#include <ThunderEgg/Vector.h>
#include <valarray>
namespace ThunderEgg
{
/**
 * @brief Vector class the uses std::valarray for data storage
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class ValVector : public Vector<D>

{
	private:
	/**
	 * @brief the number of cells in each patch
	 */
	int patch_stride;
	/**
	 * @brief the number of non-ghost cells in each direction of the patch
	 */
	std::array<int, D> lengths;
	/**
	 * @brief the strides for each axis of the patch
	 */
	std::array<int, D> strides;
	/**
	 * @brief the number of ghost cells on each side of the patch
	 */
	int num_ghost_cells;
	/**
	 * @brief the offset of the first element in each patch
	 */
	int first_offset;
	/**
	 * @brief the underlying vector
	 */
	std::valarray<double> vec;

	/**
	 * @brief Calculate the number of local (non-ghost) cells
	 *
	 * @param lengths the number of cells in each direction
	 * @param num_patches the number of patches in this vector
	 * @return int the number of local (non-ghost) cells
	 */
	static int GetNumLocalCells(const std::array<int, D> &lengths, int num_patches)
	{
		int num_cells_in_patch = 1;
		for (size_t i = 0; i < D; i++) {
			num_cells_in_patch *= lengths[i];
		}
		return num_patches * num_cells_in_patch;
	}

	public:
	/**
	 * @brief Construct a new ValVector object
	 *
	 * @param comm the MPI_Comm to use
	 * @param lengths the nubmer of (non-ghost) cells in each direction of a patch
	 * @param num_ghost_cells the number of ghost cells padding each side of a patch
	 * @param num_components the number of components for each cell
	 * @param num_patches the number of patches in this vector
	 */
	ValVector(MPI_Comm comm, const std::array<int, D> &lengths, int num_ghost_cells,
	          int num_components, int num_patches)
	: Vector<D>(comm, num_components, num_patches, GetNumLocalCells(lengths, num_patches)),
	  lengths(lengths), num_ghost_cells(num_ghost_cells)
	{
		int size            = num_components;
		int my_first_offset = 0;
		for (size_t i = 0; i < D; i++) {
			strides[i] = size;
			size *= (this->lengths[i] + 2 * num_ghost_cells);
			my_first_offset += strides[i] * num_ghost_cells;
		}
		first_offset = my_first_offset;
		patch_stride = size;
		size *= num_patches;
		vec.resize(size);
	}
	/**
	 * @brief Get a new ValVector object for a given Domain
	 *
	 * @param domain the Domain
	 * @param num_components the number of components for each cell
	 * @return std::shared_ptr<ValVector<D>> the new Vector
	 */
	static std::shared_ptr<ValVector<D>> GetNewVector(std::shared_ptr<const Domain<D>> domain,
	                                                  int num_components)
	{
		return std::shared_ptr<ValVector<D>>(
		new ValVector<D>(MPI_COMM_WORLD, domain->getNs(), domain->getNumGhostCells(),
		                 num_components, domain->getNumLocalPatches()));
	}
	LocalData<D> getLocalData(int component_index, int local_patch_index) override
	{
		double *data = &vec[patch_stride * local_patch_index + first_offset + component_index];
		return LocalData<D>(data, strides, lengths, num_ghost_cells, nullptr);
	}
	const LocalData<D> getLocalData(int component_index, int local_patch_index) const override
	{
		double *data = const_cast<double *>(
		&vec[patch_stride * local_patch_index + component_index + first_offset + component_index]);
		return LocalData<D>(data, strides, lengths, num_ghost_cells, nullptr);
	}

	/**
	 * @brief Get the number of ghost cells padding each side of the patches
	 *
	 * @return int the number of ghost cells padding each side of the patches
	 */
	int getNumGhostCells() const
	{
		return num_ghost_cells;
	}
	/**
	 * @brief Get a reference to the underlying valarray
	 *
	 * @return std::valarray<int>& a reference to the underlying valarray
	 */
	std::valarray<double> &getValArray()
	{
		return vec;
	}
};
} // namespace ThunderEgg
#endif
