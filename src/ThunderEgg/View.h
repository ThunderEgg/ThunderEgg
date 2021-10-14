/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#ifndef THUNDEREGG_VIEW_H
#define THUNDEREGG_VIEW_H
/**
 * @file
 *
 * @brief View class
 */
#include <ThunderEgg/Config.h>
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/RuntimeError.h>
#include <array>
#include <memory>
namespace ThunderEgg
{
/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <typename T, int D> class View
{
	public:
	using T_ptr = typename std::add_pointer<T>::type;

	private:
	/**
	 * @brief Pointer to the first non-ghost cell value
	 */
	T_ptr data = nullptr;
	/**
	 * @brief the strides between each element, in each direction
	 */
	std::array<int, D> strides;
	/**
	 * @brief the coordianate of the first ghost element
	 */
	std::array<int, D> ghost_start;
	/**
	 * @brief the coordinate of the first non-ghost element
	 */
	std::array<int, D> start;
	/**
	 * @brief the coordinate of the last non-ghost
	 */
	std::array<int, D> end;
	/**
	 * @brief the corrdinate of the last ghost element
	 */
	std::array<int, D> ghost_end;

	/**
	 * @brief check that a coordinate is in bounds
	 * will throw an exception if it isn't
	 *
	 * @param coord the coordinate
	 */
	void checkCoordIsInBounds(const std::array<int, D> &coord) const
	{
		bool is_valid_coord = true;
		for (int i = 0; i < D; i++) {
			is_valid_coord = is_valid_coord && (coord[i] >= ghost_start[i] && coord[i] <= ghost_end[i]);
		}
		if (!is_valid_coord) {
			// oob coord
			throw RuntimeError("index for view is out of bounds");
		}
	}

	protected:
	inline int getIndex(const std::array<int, D> &coord) const
	{
		int idx = 0;
		Loop::Unroll<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
		return idx;
	}

	T_ptr getData() const
	{
		return data;
	}

	public:
	/**
	 * @brief Constructs a view of size 0
	 */
	View()
	{
		strides.fill(0);
		ghost_start.fill(0);
		start.fill(0);
		end.fill(-1);
		ghost_end.fill(-1);
	}
	/**
	 * @brief Construct a new View object
	 *
	 * @param data pointer to the first element in the patch (non-ghost cell element)
	 * @param strides the strides in each direction
	 * @param lengths the lengths in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @param ldm the local data manager for the data
	 */
	View(T_ptr                     data,
	     const std::array<int, D> &strides,
	     const std::array<int, D> &ghost_start,
	     const std::array<int, D> &start,
	     const std::array<int, D> &end,
	     const std::array<int, D> &ghost_end)
	: strides(strides),
	  ghost_start(ghost_start),
	  start(start),
	  end(end),
	  ghost_end(ghost_end)
	{
		this->data = data - this->getIndex(ghost_start);
	}

	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * @param coord the coordinate
	 * @return double& the element
	 */
	inline T &operator[](const std::array<int, D> &coord) const
	{
		if constexpr (ENABLE_DEBUG) {
			this->checkCoordIsInBounds(coord);
		}
		return data[this->getIndex(coord)];
	}

	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * for example uasage in 3d is view(ix,iy,iz)
	 *
	 * @tparam Types the types
	 * @param args the coordnate in x,y,z form
	 * @return const double&  the value
	 */
	template <class... Types> inline T &operator()(Types... args) const
	{
		static_assert(sizeof...(args) == D, "incorrect number of arguments");
		if constexpr (ENABLE_DEBUG) {
			this->checkCoordIsInBounds({args...});
		}
		return data[this->getIndex({args...})];
	}

	/**
	 * @brief Set the value at a coordinate to the specified value
	 *
	 * @param coord the coordinate
	 * @param value the value
	 */
	inline void set(const std::array<int, D> &coord, T value) const
	{
		if constexpr (std::is_const<T>::value) {
			if constexpr (ENABLE_DEBUG) {
				// check that only ghost cells are being modified
				bool is_interior_coord = true;
				for (int i = 0; i < D; i++) {
					is_interior_coord = is_interior_coord && (coord[i] >= start[i] && coord[i] <= end[i]);
				}
				if (is_interior_coord) {
					// intertior coord
					throw RuntimeError("interior value of const view is being modified");
				}
				checkCoordIsInBounds(coord);
			}
			int idx = 0;
			Loop::Unroll<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
			const_cast<double *>(data)[idx] = value;
		} else {
			if constexpr (ENABLE_DEBUG) {
				this->checkCoordIsInBounds(coord);
			}
			data[this->getIndex(coord)] = value;
		}
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

	operator View<std::add_const_t<T>, D>() const
	{
		return View<std::add_const_t<T>, D>(data + getIndex(getGhostStart()), strides, ghost_start, start, end, ghost_end);
	}
};
extern template class View<double, 1>;
extern template class View<double, 2>;
extern template class View<double, 3>;
extern template class View<double, 4>;
extern template class View<const double, 1>;
extern template class View<const double, 2>;
extern template class View<const double, 3>;
extern template class View<const double, 4>;
} // namespace ThunderEgg
#endif