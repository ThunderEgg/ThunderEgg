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

#ifndef THUNDEREGG_VIEW_H
#define THUNDEREGG_VIEW_H
#include <ThunderEgg/ConstView.h>
namespace ThunderEgg
{
/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <int D> class View : public ConstView<D>
{
	private:
	/**
	 * @brief Pointer to the first non-ghost cell value
	 */
	double *data = nullptr;

	public:
	/**
	 * @brief Constructs a view of size 0
	 */
	View() : ConstView<D>() {}
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
	: ConstView<D>(data, strides, ghost_start, start, end, ghost_end, ldm)
	{
		this->data = data - this->getIndex(ghost_start);
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
	template <class... Types> inline double &operator()(Types... args)
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
	inline void set(const std::array<int, D> &coord, double value)
	{
		if constexpr (ENABLE_DEBUG) {
			this->checkCoordIsInBounds(coord);
		}
		data[this->getIndex(coord)] = value;
	}

	using ConstView<D>::operator[];
	using ConstView<D>::operator();
	using ConstView<D>::set;
};
extern template class View<1>;
extern template class View<2>;
extern template class View<3>;
extern template class View<4>;
} // namespace ThunderEgg
#endif