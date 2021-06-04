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

#ifndef THUNDEREGG_PATCHARRAY_H
#define THUNDEREGG_PATCHARRAY_H
#include <ThunderEgg/View.h>
namespace ThunderEgg
{
/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <int D> class PatchArray
{
	private:
	/**
	 * @brief data
	 */
	std::vector<double> data;
	/**
	 * @brief the strides between each element, in each direction
	 */
	std::array<int, D> strides;
	/**
	 * @brief The number of cells in each direction
	 */
	std::array<int, D> lengths;
	/**
	 * @brief the coordinate of the first non-ghost element
	 */
	std::array<int, D> start;
	/**
	 * @brief the coordinate of the last non-ghost
	 */
	std::array<int, D> end;
	/**
	 * @brief the coordianate of the first ghost element
	 */
	std::array<int, D> ghost_start;
	/**
	 * @brief the corrdinate of the last ghost element
	 */
	std::array<int, D> ghost_end;
	/**
	 * @brief The ViewManager that does any necessary cleanup
	 */
	std::shared_ptr<ViewManager> ldm;
	/**
	 * @brief The number of ghost cells on each side of the patch
	 */
	int num_ghost_cells;
	/**
	 * @brief Starting index of first non-ghost cell
	 */
	int start_idx;

	/**
	 * @brief Get a slice for a slide
	 *
	 * @param s the side of patch for the slice
	 * @param offset the offset, with {0, 0} being the first slice of non-ghost cell values, and {-1, -1} being
	 * the first slice of ghost cell values
	 * @return View<1> the resulting slice
	 */
	template <int M> View<M> getSliceOnPriv(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		std::array<int, M>         new_strides;
		std::array<int, M>         new_lengths;
		double *                   new_data      = const_cast<double *>(data.data()) + start_idx;
		std::array<Side<D>, D - M> sides         = f.getSides();
		size_t                     lengths_index = 0;
		size_t                     sides_index   = 0;
		for (size_t axis = 0; axis < (size_t) D; axis++) {
			if (sides_index < sides.size() && sides[sides_index].getAxisIndex() == axis) {
				if (sides[sides_index].isLowerOnAxis()) {
					new_data += offset[sides_index] * strides[axis];
				} else {
					new_data += (lengths[axis] - 1 - offset[sides_index]) * strides[axis];
				}
				sides_index++;
			} else {
				new_lengths[lengths_index] = lengths[axis];
				new_strides[lengths_index] = strides[axis];
				lengths_index++;
			}
		}
		return View<M>(new_data, new_strides, new_lengths, num_ghost_cells, ldm);
	}

	public:
	/**
	 * @brief Construct a new View object
	 *
	 * @param lengths the lengths in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 */
	PatchArray(const std::array<int, D> &lengths, int num_ghost_cells) : lengths(lengths), num_ghost_cells(num_ghost_cells)
	{
		start.fill(0);
		end = lengths;
		for (int i = 0; i < D; i++) {
			end[i]--;
		}
		ghost_start = start;
		ghost_end   = end;
		for (size_t i = 0; i < D; i++) {
			ghost_start[i] -= num_ghost_cells;
			ghost_end[i] += num_ghost_cells;
		}
		int size  = 1;
		start_idx = 0;
		for (int i = 0; i < D; i++) {
			strides[i] = size;
			start_idx += num_ghost_cells * strides[i];
			size *= lengths[i] + 2 * num_ghost_cells;
		}
		data.resize(size);
	}
	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * @param coord the coordinate
	 * @return double& the element
	 */
	inline double &operator[](const std::array<int, D> &coord)
	{
		int idx = start_idx;
		loop<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
		return data[idx];
	}
	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * @param coord the coordinate
	 * @return double& the element
	 */
	inline const double &operator[](const std::array<int, D> &coord) const
	{
		int idx = start_idx;
		loop<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
		return data[idx];
	}
	/**
	 * @brief Get a slice in a given face
	 *
	 * @param f the face
	 * @param offset the offset of the value {0,..,0} is the non ghost cell touching the corner. {-1,..,-1} is the first ghost cell touching that
	 * face
	 * @return double& the value
	 */
	template <int M> View<M> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset)
	{
		return getSliceOnPriv(f, offset);
	}
	/**
	 * @brief Get a slice in a given face
	 *
	 * @param f the face
	 * @param offset the offset of the value {0,..,0} is the non ghost cell touching the corner. {-1,..,-1} is the first ghost cell touching that
	 * face
	 * @return double& the value
	 */
	template <int M> const View<M> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		return getSliceOnPriv(f, offset);
	}
	template <int M> View<M> getGhostSliceOn(Face<D, M> f, const std::array<size_t, D - M> &offset) const
	{
		std::array<int, M>         new_strides;
		std::array<int, M>         new_lengths;
		double *                   new_data      = const_cast<double *>(data.data()) + start_idx;
		std::array<Side<D>, D - M> sides         = f.getSides();
		size_t                     lengths_index = 0;
		size_t                     sides_index   = 0;
		for (size_t axis = 0; axis < (size_t) D; axis++) {
			if (sides_index < sides.size() && sides[sides_index].getAxisIndex() == axis) {
				if (sides[sides_index].isLowerOnAxis()) {
					new_data += (-1 - offset[sides_index]) * strides[axis];
				} else {
					new_data += (lengths[axis] + offset[sides_index]) * strides[axis];
				}
				sides_index++;
			} else {
				new_lengths[lengths_index] = lengths[axis];
				new_strides[lengths_index] = strides[axis];
				lengths_index++;
			}
		}
		return View<M>(new_data, new_strides, new_lengths, num_ghost_cells, ldm);
	}
	/**
	 * @brief Get the Lengths of the patch in each direction
	 */
	const std::array<int, D> &getLengths() const
	{
		return lengths;
	}
	/**
	 * @brief Get the strides of the patch in each direction
	 */
	const std::array<int, D> &getStrides() const
	{
		return strides;
	}
	/**
	 * @brief Get the coordinate of the first element
	 */
	const std::array<int, D> &getStart() const
	{
		return start;
	}
	/**
	 * @brief Get the coordinate of the last element
	 */
	const std::array<int, D> &getEnd() const
	{
		return end;
	}
	/**
	 * @brief Get the coordinate of the first ghost cell element
	 */
	const std::array<int, D> &getGhostStart() const
	{
		return ghost_start;
	}
	/**
	 * @brief Get the coordinate of the last ghost cell element
	 */
	const std::array<int, D> &getGhostEnd() const
	{
		return ghost_end;
	}
	/**
	 * @brief Get the number of ghost cells on each side of the patch
	 */
	int getNumGhostCells() const
	{
		return num_ghost_cells;
	}
};
extern template class PatchArray<1>;
extern template class PatchArray<2>;
extern template class PatchArray<3>;
} // namespace ThunderEgg
#endif