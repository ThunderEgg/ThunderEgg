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

#include <ThunderEgg/Vector.h>
namespace ThunderEgg
{
/**
 * @brief Check if a coordinate is a ghost coordinate
 *
 * @tparam D the number of Cartesian Dimensions
 * @param coord the coordinate
 * @param ns the number of cells in each direction of the patch
 * @param num_ghost_cells the number of ghost cells on each side of the patch
 * @return true if it is a ghost coordinate
 * @return false if it is not a ghost coordinate
 */
template <size_t D>
bool isGhost(const std::array<int, D> &coord, const std::array<int, D> &ns, int num_ghost_cells)
{
	for (size_t i = 0; i < D; i++) {
		if (coord[i] < 0 || coord[i] >= ns[i]) {
			return true;
		}
	}
	return false;
}
template <int D>
class MockVector : public Vector<D>
{
	public:
	int GetNumLocalCells(int num_local_patches, std::array<int, D> ns)
	{
		int num_cells = 1;
		for (size_t i = 0; i < D; i++) {
			num_cells *= ns[i];
		}
		return num_cells * num_local_patches;
	}
	/**
	 * @brief underlying data storage
	 */
	std::vector<double> data;
	/**
	 * @brief vector of View objects
	 */
	std::vector<ComponentView<D>> local_data;
	/**
	 * @brief Construct a new MockVector object
	 *
	 * @param comm the MPI_Comm
	 * @param num_components number of components for each cell
	 * @param num_local_patches number of local patches
	 * @param num_ghost_cells number of ghost cells
	 * @param ns the number of cells in each direction of a patch
	 */
	MockVector(MPI_Comm comm, int num_components, int num_local_patches, int num_ghost_cells,
	           std::array<int, D> ns)
	: Vector<D>(comm, num_components, num_local_patches, GetNumLocalCells(num_local_patches, ns))
	{
		std::array<int, 3> strides;
		int                patch_stride = num_components;
		int                first_offset = 0;
		for (size_t i = 0; i < 3; i++) {
			strides[i] = patch_stride;
			patch_stride *= (ns[i] + 2 * num_ghost_cells);
			first_offset += strides[i] * num_ghost_cells;
		}
		int size = patch_stride * num_local_patches;

		data.resize(size);

		local_data.reserve(num_local_patches);

		for (int i = 0; i < num_local_patches; i++) {
			for (int c = 0; c < num_components; c++) {
				double *data_ptr = data.data() + i * patch_stride + first_offset + c;
				local_data.emplace_back(data_ptr, strides, ns, num_ghost_cells);
			}
		}
	}
	ComponentView<D> getComponentView(int component_index, int patch_local_index) override
	{
		return local_data[patch_local_index * this->getNumComponents() + component_index];
	}
	const ComponentView<D> getComponentView(int component_index, int patch_local_index) const override
	{
		return local_data[patch_local_index * this->getNumComponents() + component_index];
	}
};
} // namespace ThunderEgg