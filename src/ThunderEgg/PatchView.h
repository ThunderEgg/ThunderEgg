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
template <typename T, int D> class PatchView : public View<T, D + 1>
{
	private:
	template <int M> struct SliceInfo {
		std::array<int, M> strides;
		std::array<int, M> ghost_start;
		std::array<int, M> start;
		std::array<int, M> end;
		std::array<int, M> ghost_end;
		std::array<int, D> first_value;
	};

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
	template <int M> SliceInfo<M> getSliceOnPriv(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		SliceInfo<M> info;

		info.first_value = this->getGhostStart();

		std::array<Side<D>, D - M> sides         = f.getSides();
		size_t                     lengths_index = 0;
		size_t                     sides_index   = 0;

		for (size_t axis = 0; axis < (size_t) D; axis++) {
			if (sides_index < sides.size() && sides[sides_index].getAxisIndex() == axis) {
				if (sides[sides_index].isLowerOnAxis()) {
					info.first_value[axis] = offset[sides_index];
				} else {
					info.first_value[axis] = this->getEnd()[axis] - offset[sides_index];
				}
				sides_index++;
			} else {
				info.strides[lengths_index]     = this->getStrides()[axis];
				info.ghost_start[lengths_index] = this->getGhostStart()[axis];
				info.start[lengths_index]       = this->getStart()[axis];
				info.end[lengths_index]         = this->getEnd()[axis];
				info.ghost_end[lengths_index]   = this->getGhostEnd()[axis];
				lengths_index++;
			}
		}
		return info;
	}

	public:
	using T_ptr = typename View<T, D + 1>::T_ptr;
	/**
	 * @brief Construct a new PatchView with size 0
	 */
	PatchView() : View<T, D + 1>()
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
	PatchView(T_ptr                              data,
	          const std::array<int, D + 1> &     strides,
	          const std::array<int, D + 1> &     lengths,
	          int                                num_ghost_cells,
	          std::shared_ptr<const ViewManager> ldm = nullptr)
	: View<T, D + 1>(data,
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
	 * @brief Get the slice on a given face
	 *
	 * @tparam M the dimension of the face
	 * @param f the face
	 * @param offset offset the offset of the value {0,..,0} is the non ghost cell touching the face. {-1,..,-1} is the first ghost cell touching that
	 * face
	 * @return View<M> a view to the slice on the face
	 */
	template <int M> View<T, M> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		SliceInfo<M> info = getSliceOnPriv<M>(f, offset);

		T_ptr new_data = (&(*this)[info.first_value]);
		return View<T, M>(new_data, info.strides, info.ghost_start, info.start, info.end, info.ghost_end, this->getPatchViewDataManager());
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
	template <int M> View<typename std::remove_const<T>::type, M> getGhostSliceOn(Face<D, M> f, const std::array<size_t, D - M> &offset) const
	{
		using noconst_T     = typename std::remove_const<T>::type;
		using noconst_T_ptr = typename std::add_pointer<noconst_T>::type;
		std::array<int, M> new_strides;
		std::array<int, M> new_ghost_start;
		std::array<int, M> new_start;
		std::array<int, M> new_end;
		std::array<int, M> new_ghost_end;
		std::array<int, D> first_value;

		first_value = this->getGhostStart();

		std::array<Side<D>, D - M> sides         = f.getSides();
		size_t                     lengths_index = 0;
		size_t                     sides_index   = 0;

		for (size_t axis = 0; axis < (size_t) D; axis++) {
			if (sides_index < sides.size() && sides[sides_index].getAxisIndex() == axis) {
				if (sides[sides_index].isLowerOnAxis()) {
					first_value[axis] = -1 - offset[sides_index];
				} else {
					first_value[axis] = this->getEnd()[axis] + 1 + offset[sides_index];
				}
				sides_index++;
			} else {
				new_strides[lengths_index]     = this->getStrides()[axis];
				new_ghost_start[lengths_index] = this->getGhostStart()[axis];
				new_start[lengths_index]       = this->getStart()[axis];
				new_end[lengths_index]         = this->getEnd()[axis];
				new_ghost_end[lengths_index]   = this->getGhostEnd()[axis];
				lengths_index++;
			}
		}

		noconst_T_ptr new_data = const_cast<noconst_T_ptr>(&(*this)[first_value]); // Thunderegg doesn't care if values in ghosts are modified
		return View<noconst_T, M>(new_data, new_strides, new_ghost_start, new_start, new_end, new_ghost_end, this->getPatchViewDataManager());
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
} // namespace ThunderEgg
#endif