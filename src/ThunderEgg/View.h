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

#ifndef THUNDEREGG_MULTIDIMENSIONALVIEW_H
#define THUNDEREGG_MULTIDIMENSIONALVIEW_H
#include <ThunderEgg/Config.h>
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/ViewManager.h>
#include <array>
#include <memory>
namespace ThunderEgg
{
/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <int D> class View
{
	private:
	/**
	 * @brief Pointer to the first non-ghost cell value
	 */
	double *data = nullptr;
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
	 * @brief The ViewManager that does any necessary cleanup
	 */
	std::shared_ptr<const ViewManager> ldm;

	template <int idx, class Type, class... Types> int getIndex(Type t, Types... args) const
	{
		return strides[idx] * t + getIndex<idx + 1>(args...);
	}
	template <int idx, class Type> int getIndex(Type t) const
	{
		return strides[idx] * t;
	}

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
	View(double *                           data,
	     const std::array<int, D> &         strides,
	     const std::array<int, D> &         ghost_start,
	     const std::array<int, D> &         start,
	     const std::array<int, D> &         end,
	     const std::array<int, D> &         ghost_end,
	     std::shared_ptr<const ViewManager> ldm = nullptr)
	: data(data),
	  strides(strides),
	  ghost_start(ghost_start),
	  start(start),
	  end(end),
	  ghost_end(ghost_end),
	  ldm(ldm)
	{
	}

	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * @param coord the coordinate
	 * @return double& the element
	 */
	inline double &operator[](const std::array<int, D> &coord)
	{
		if constexpr (ENABLE_DEBUG) {
			checkCoordIsInBounds(coord);
		}
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
		if constexpr (ENABLE_DEBUG) {
			checkCoordIsInBounds(coord);
		}
		int idx = 0;
		loop<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
		return data[idx];
	}
	template <class... Types> inline double &operator()(Types... args)
	{
		static_assert(sizeof...(args) == D, "incorrect number of arguments");
		if constexpr (ENABLE_DEBUG) {
			checkCoordIsInBounds({args...});
		}
		return data[getIndex<0>(args...)];
	}
	template <class... Types> inline const double &operator()(Types... args) const
	{
		static_assert(sizeof...(args) == D, "incorrect number of arguments");
		if constexpr (ENABLE_DEBUG) {
			checkCoordIsInBounds({args...});
		}
		return data[getIndex<0>(args...)];
	}
	inline void set(const std::array<int, D> &coord, double value)
	{
		if constexpr (ENABLE_DEBUG) {
			checkCoordIsInBounds(coord);
		}
		int idx = 0;
		loop<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
		data[idx] = value;
	}
	inline void set(const std::array<int, D> &coord, double value) const
	{
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
		loop<0, D - 1>([&](int i) { idx += strides[i] * coord[i]; });
		data[idx] = value;
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
	const std::shared_ptr<const ViewManager> getComponentViewDataManager() const
	{
		return ldm;
	}
};
extern template class View<1>;
extern template class View<2>;
extern template class View<3>;
extern template class View<4>;
} // namespace ThunderEgg
#endif