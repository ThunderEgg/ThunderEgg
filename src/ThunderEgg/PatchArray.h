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
/**
 * @file
 *
 * @brief PatchArray class
 */
#include <ThunderEgg/PatchView.h>
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
	std::vector<double>  vector;
	PatchView<double, D> view;

	public:
	/**
	 * @brief Construct a new PatchArray object of size zero
	 *
	 */
	PatchArray() = default;
	/**
	 * @brief Construct a new PatchArray object
	 *
	 * @param lengths the lengths in each direction
	 * @param num_components the number of components
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 */
	PatchArray(const std::array<int, D> &lengths, int num_components, int num_ghost_cells)
	{
		std::array<int, D + 1> view_lengths;
		for (int i = 0; i < D; i++) {
			view_lengths[i] = lengths[i];
		}
		view_lengths[D] = num_components;
		std::array<int, D + 1> strides;
		strides[0] = 1;
		for (int i = 1; i < D + 1; i++) {
			strides[i] = strides[i - 1] * (lengths[i - 1] + 2 * num_ghost_cells);
		}
		int size = strides[D] * num_components;
		vector.resize(size);
		view = PatchView<double, D>(vector.data(), strides, view_lengths, num_ghost_cells);
	}

	/**
	 * @brief Copy constructor
	 *
	 * @param other the array to copy
	 */
	PatchArray(const PatchArray<D> &other)
	: vector(other.vector),
	  view(vector.data(), other.getStrides(), other.getGhostStart(), other.getStart(), other.getEnd(), other.getGhostEnd())

	{
	}

	/**
	 * @brief Copy assignment
	 *
	 * @param other the array to copy
	 * @return PatchArray<D>&  this
	 */
	PatchArray<D> &operator=(const PatchArray<D> &other)
	{
		vector = other.vector;
		view = PatchView<double, D>(vector.data(), other.getStrides(), other.getGhostStart(), other.getStart(), other.getEnd(), other.getGhostEnd());
		return *this;
	}

	/**
	 * @brief Get the slice on a given face
	 *
	 * @tparam M the dimension of the face
	 * @param f the face
	 * @param offset offset the offset of the value {0,..,0} is the non ghost cell touching the face. {-1,..,-1} is the first ghost cell touching that
	 * face
	 * @return View<M+1> a view to the slice on the face
	 */
	template <int M> inline View<double, M + 1> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset)
	{
		return view.template getSliceOn<M>(f, offset);
	}

	/**
	 * @brief Get the slice on a given face
	 *
	 * @tparam M the dimension of the face
	 * @param f the face
	 * @param offset offset the offset of the value {0,..,0} is the non ghost cell touching the face. {-1,..,-1} is the first ghost cell touching that
	 * face
	 * @return ConstView<M> a view to the slice on the face
	 */
	template <int M> inline View<const double, M + 1> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		return View<const double, M + 1>(view.template getSliceOn<M>(f, offset));
	}

	/**
	 * @brief Get the gosts slice on a given face
	 *
	 * @tparam M the dimension of the face
	 * @param f the face
	 * @param offset offset the offset of the value {0,..,0} first ghost cell slice touching that face
	 * face
	 * @return View<M> a view to the slice on the face
	 */
	template <int M> inline View<double, M + 1> getGhostSliceOn(Face<D, M> f, const std::array<size_t, D - M> &offset) const
	{
		return view.template getGhostSliceOn<M>(f, offset);
	}

	inline const double &operator[](const std::array<int, D + 1> &coord) const
	{
		return view[coord];
	}
	template <class... Types> inline const double &operator()(Types... args) const
	{
		return view(args...);
	}
	inline void set(const std::array<int, D + 1> &coord, double value) const
	{
		view.set(coord, value);
	}
	inline double &operator[](const std::array<int, D + 1> &coord)
	{
		return view[coord];
	}
	template <class... Types> inline double &operator()(Types... args)
	{
		return view(args...);
	}
	inline void set(const std::array<int, D + 1> &coord, double value)
	{
		view.set(coord, value);
	}

	/**
	 * @brief Get the strides of the patch in each direction
	 */
	inline const std::array<int, D + 1> &getStrides() const
	{
		return view.getStrides();
	}
	/**
	 * @brief Get the coordinate of the first element
	 */
	inline const std::array<int, D + 1> &getStart() const
	{
		return view.getStart();
	}
	/**
	 * @brief Get the coordinate of the last element
	 */
	inline const std::array<int, D + 1> &getEnd() const
	{
		return view.getEnd();
	}
	/**
	 * @brief Get the coordinate of the first ghost cell element
	 */
	inline const std::array<int, D + 1> &getGhostStart() const
	{
		return view.getGhostStart();
	}
	/**
	 * @brief Get the coordinate of the last ghost cell element
	 */
	inline const std::array<int, D + 1> &getGhostEnd() const
	{
		return view.getGhostEnd();
	}
	/**
	 * @brief Get the View for the array
	 *
	 * @return const PatchView<double, D>& the View
	 */
	const PatchView<double, D> &getView()
	{
		return view;
	}
};
extern template class PatchArray<1>;
extern template class PatchArray<2>;
extern template class PatchArray<3>;
} // namespace ThunderEgg
#endif