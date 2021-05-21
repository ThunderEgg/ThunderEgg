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

#ifndef THUNDEREGG_LOCALDATA_H
#define THUNDEREGG_LOCALDATA_H
#include <ThunderEgg/Face.h>
#include <ThunderEgg/LocalDataManager.h>
#include <ThunderEgg/Loops.h>
#include <memory>
namespace ThunderEgg
{
/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <int D> class LocalData
{
	private:
	/**
	 * @brief Pointer to the first non-ghost cell value
	 */
	double *data;
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
	 * @brief The LocalDataManager that does any necessary cleanup
	 */
	std::shared_ptr<LocalDataManager> ldm;
	/**
	 * @brief The number of ghost cells on each side of the patch
	 */
	int num_ghost_cells;

	/**
	 * @brief Get a slice for a slide
	 *
	 * @param s the side of patch for the slice
	 * @param offset the offset, with {0, 0} being the first slice of non-ghost cell values, and {-1, -1} being
	 * the first slice of ghost cell values
	 * @return LocalData<1> the resulting slice
	 */
	template <int M> LocalData<M> getSliceOnPriv(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		std::array<int, M>         new_strides;
		std::array<int, M>         new_lengths;
		double *                   new_data      = data;
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
		return LocalData<M>(new_data, new_strides, new_lengths, 0, ldm);
	}

	template <int idx, class Type, class... Types> int getIndex(Type t, Types... args) const
	{
		return strides[idx] * t + getIndex<idx + 1>(args...);
	}
	template <int idx, class Type> int getIndex(Type t) const
	{
		return strides[idx] * t;
	}

	public:
	LocalData() {}
	/**
	 * @brief Construct a new LocalData object
	 *
	 * @param data_in pointer to the first element in the patch (non-ghost cell element)
	 * @param strides_in the strides in each direction
	 * @param lengths_in the lengths in each direction
	 * @param num_ghost_cells_in the number of ghost cells on each side of the patch
	 * @param ldm_in the local data manager for the data
	 */
	LocalData(double *                          data_in,
	          const std::array<int, D> &        strides_in,
	          const std::array<int, D> &        lengths_in,
	          int                               num_ghost_cells_in,
	          std::shared_ptr<LocalDataManager> ldm_in = nullptr)
	: data(data_in),
	  strides(strides_in),
	  lengths(lengths_in),
	  ldm(ldm_in),
	  num_ghost_cells(num_ghost_cells_in)
	{
		start.fill(0);
		end = lengths;
		for (size_t i = 0; i < D; i++) {
			start[i] = 0;
			end[i]--;
		}
		ghost_start = start;
		ghost_end   = end;
		for (size_t i = 0; i < D; i++) {
			ghost_start[i] -= num_ghost_cells;
			ghost_end[i] += num_ghost_cells;
		}
	}
	/**
	 * @brief Get the pointer the data at the specified coordinate
	 *
	 * @param coord the coordianate
	 * @return double* the pointer
	 */
	inline double *getPtr(const std::array<int, D> &coord)
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data + idx;
	}
	/**
	 * @brief Get the pointer the data at the specified coordinate
	 *
	 * @param coord the coordianate
	 * @return double* the pointer
	 */
	inline const double *getPtr(const std::array<int, D> &coord) const
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data + idx;
	}
	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * @param coord the coordinate
	 * @return double& the element
	 */
	inline double &operator[](const std::array<int, D> &coord)
	{
		int idx = 0;
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
		int idx = 0;
		loop<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
		return data[idx];
	}
	template <class... Types> inline double &operator()(Types... args)
	{
		static_assert(sizeof...(args) == D, "incorrect number of arguments");
		return data[getIndex<0>(args...)];
	}
	template <class... Types> inline const double &operator()(Types... args) const
	{
		static_assert(sizeof...(args) == D, "incorrect number of arguments");
		return data[getIndex<0>(args...)];
	}
	/**
	 * @brief Get a slice in a given face
	 *
	 * @param f the face
	 * @param offset the offset of the value {0,..,0} is the non ghost cell touching the corner. {-1,..,-1} is the first ghost cell touching that
	 * face
	 * @return double& the value
	 */
	template <int M> LocalData<M> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset)
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
	template <int M> const LocalData<M> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		return getSliceOnPriv(f, offset);
	}
	template <int M> LocalData<M> getGhostSliceOn(Face<D, M> f, const std::array<size_t, D - M> &offset) const
	{
		std::array<int, M>         new_strides;
		std::array<int, M>         new_lengths;
		double *                   new_data      = data;
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
		return LocalData<M>(new_data, new_strides, new_lengths, 0, ldm);
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
	/**
	 * @brief Get the pointer to the first element
	 */
	double *getPtr() const
	{
		return data;
	}
};
extern template class LocalData<1>;
extern template class LocalData<2>;
extern template class LocalData<3>;
/*
template <> inline const double &LocalData<2>::operator[](const std::array<int, 2> &coord) const
{
    return data[strides[0] * coord[0] + strides[1] * coord[1]];
}
template <> inline double &LocalData<2>::operator[](const std::array<int, 2> &coord)
{
    return data[strides[0] * coord[0] + strides[1] * coord[1]];
}
*/
} // namespace ThunderEgg
#endif