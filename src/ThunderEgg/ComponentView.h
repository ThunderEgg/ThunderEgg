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

#ifndef THUNDEREGG_COMPONENTVIEW_H
#define THUNDEREGG_COMPONENTVIEW_H
#include <ThunderEgg/Config.h>
#include <ThunderEgg/Face.h>
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/View.h>
#include <ThunderEgg/ViewManager.h>
#include <memory>
namespace ThunderEgg
{
/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <int D> class ComponentView : public View<D>
{
	private:
	/**
	 * @brief The number of cells in each direction
	 */
	std::array<int, D> lengths;
	/**
	 * @brief The number of ghost cells on each side of the patch
	 */
	int num_ghost_cells = 0;

	/**
	 * @brief Return array filled with -num_ghost_cells
	 */
	template <int M> static std::array<int, M> DetermineGhostStart(int num_ghost_cells)
	{
		std::array<int, M> start;
		start.fill(-num_ghost_cells);
		return start;
	}

	/**
	 * @brief Return array filled with 0
	 */
	template <int M> static std::array<int, M> DetermineStart()
	{
		std::array<int, M> start;
		start.fill(0);
		return start;
	}

	/**
	 * @brief Return array filled with lengths-1
	 */
	template <int M> static std::array<int, M> DetermineEnd(const std::array<int, M> &lengths)
	{
		std::array<int, M> end = lengths;
		loop<0, M - 1>([&](int i) { end[i]--; });
		return end;
	}

	/**
	 * @brief Return array filled with lengths-1+num_ghost_cells
	 */
	template <int M> static std::array<int, M> DetermineGhostEnd(const std::array<int, M> &lengths, int num_ghost_cells)
	{
		std::array<int, M> end = lengths;
		loop<0, M - 1>([&](int i) { end[i] += num_ghost_cells - 1; });
		return end;
	}

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
		const std::array<int, D> &strides = this->getStrides();

		std::array<int, M> new_strides;
		std::array<int, M> new_lengths;

		std::array<int, D> first_non_ghost_value;
		first_non_ghost_value.fill(0);

		std::array<Side<D>, D - M> sides         = f.getSides();
		size_t                     lengths_index = 0;
		size_t                     sides_index   = 0;
		for (size_t axis = 0; axis < (size_t) D; axis++) {
			if (sides_index < sides.size() && sides[sides_index].getAxisIndex() == axis) {
				if (sides[sides_index].isLowerOnAxis()) {
					first_non_ghost_value[axis] = offset[sides_index];
				} else {
					first_non_ghost_value[axis] = lengths[axis] - 1 - offset[sides_index];
				}
				sides_index++;
			} else {
				new_lengths[lengths_index] = lengths[axis];
				new_strides[lengths_index] = strides[axis];
				lengths_index++;
			}
		}
		double *new_data = const_cast<double *>(&(*this)[first_non_ghost_value]);
		return View<M>(new_data,
		               new_strides,
		               DetermineGhostStart<M>(num_ghost_cells),
		               DetermineStart<M>(),
		               DetermineEnd<M>(new_lengths),
		               DetermineGhostEnd<M>(new_lengths, num_ghost_cells),
		               this->getComponentViewDataManager());
	}

	public:
	/**
	 * @brief Construct a new ComponentView with size 0
	 */
	ComponentView() : View<D>()
	{
		lengths.fill(0);
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
	ComponentView(double *                           data,
	              const std::array<int, D> &         strides,
	              const std::array<int, D> &         lengths,
	              int                                num_ghost_cells,
	              std::shared_ptr<const ViewManager> ldm = nullptr)
	: View<D>(data,
	          strides,
	          DetermineGhostStart<D>(num_ghost_cells),
	          DetermineStart<D>(),
	          DetermineEnd<D>(lengths),
	          DetermineGhostEnd<D>(lengths, num_ghost_cells),
	          ldm),
	  lengths(lengths),
	  num_ghost_cells(num_ghost_cells)
	{
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
		const std::array<int, D> &strides = this->getStrides();

		std::array<int, M> new_strides;
		std::array<int, M> new_lengths;

		std::array<int, D> first_non_ghost_value;
		first_non_ghost_value.fill(0);

		std::array<Side<D>, D - M> sides         = f.getSides();
		size_t                     lengths_index = 0;
		size_t                     sides_index   = 0;
		for (size_t axis = 0; axis < (size_t) D; axis++) {
			if (sides_index < sides.size() && sides[sides_index].getAxisIndex() == axis) {
				if (sides[sides_index].isLowerOnAxis()) {
					first_non_ghost_value[axis] = -1 - offset[sides_index];
				} else {
					first_non_ghost_value[axis] = lengths[axis] + offset[sides_index];
				}
				sides_index++;
			} else {
				new_lengths[lengths_index] = lengths[axis];
				new_strides[lengths_index] = strides[axis];
				lengths_index++;
			}
		}
		double *new_data = const_cast<double *>(&(*this)[first_non_ghost_value]);
		return View<M>(new_data,
		               new_strides,
		               DetermineGhostStart<M>(num_ghost_cells),
		               DetermineStart<M>(),
		               DetermineEnd<M>(new_lengths),
		               DetermineGhostEnd<M>(new_lengths, num_ghost_cells),
		               this->getComponentViewDataManager());
	}
	/**
	 * @brief Get the Lengths of the patch in each direction
	 */
	const std::array<int, D> &getLengths() const
	{
		return lengths;
	}
	/**
	 * @brief Get the number of ghost cells on each side of the patch
	 */
	int getNumGhostCells() const
	{
		return num_ghost_cells;
	}
};
extern template class ComponentView<1>;
extern template class ComponentView<2>;
extern template class ComponentView<3>;
/*
template <> inline const double &View<2>::operator[](const std::array<int, 2> &coord) const
{
    return data[strides[0] * coord[0] + strides[1] * coord[1]];
}
template <> inline double &View<2>::operator[](const std::array<int, 2> &coord)
{
    return data[strides[0] * coord[0] + strides[1] * coord[1]];
}
*/
} // namespace ThunderEgg
#endif