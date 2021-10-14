/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

#ifndef THUNDEREGG_PATCHVIEW_H
#define THUNDEREGG_PATCHVIEW_H
#include <ThunderEgg/ComponentView.h>
/**
 * @file
 *
 * @brief PatchView class
 */

namespace ThunderEgg
{
/**
 * @brief View for accessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <typename T, int D> class PatchView : public View<T, D + 1>
{
	private:
	template <int M> struct SliceInfo {
		std::array<int, M + 1> strides;
		std::array<int, M + 1> ghost_start;
		std::array<int, M + 1> start;
		std::array<int, M + 1> end;
		std::array<int, M + 1> ghost_end;
		std::array<int, D + 1> first_value;
	};

	/**
	 * @brief Return array filled with -num_ghost_cells
	 */
	template <int M> static std::array<int, M + 1> DetermineGhostStart(int num_ghost_cells)
	{
		std::array<int, M + 1> start;
		start.fill(-num_ghost_cells);
		start[M] = 0;
		return start;
	}

	/**
	 * @brief Return array filled with 0
	 */
	template <int M> static std::array<int, M + 1> DetermineStart()
	{
		std::array<int, M + 1> start;
		start.fill(0);
		return start;
	}

	/**
	 * @brief Return array filled with lengths-1
	 */
	template <int M> static std::array<int, M + 1> DetermineEnd(const std::array<int, M + 1> &lengths)
	{
		std::array<int, M + 1> end = lengths;
		Loop::Unroll<0, M>([&](int i) { end[i]--; });
		return end;
	}

	/**
	 * @brief Return array filled with lengths-1+num_ghost_cells
	 */
	template <int M> static std::array<int, M + 1> DetermineGhostEnd(const std::array<int, M + 1> &lengths, int num_ghost_cells)
	{
		std::array<int, M + 1> end = lengths;
		Loop::Unroll<0, M - 1>([&](int i) { end[i] += num_ghost_cells - 1; });
		end[M] = lengths[M] - 1;
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
		info.strides[M]     = this->getStrides()[D];
		info.ghost_start[M] = this->getGhostStart()[D];
		info.start[M]       = this->getStart()[D];
		info.end[M]         = this->getEnd()[D];
		info.ghost_end[M]   = this->getGhostEnd()[D];
		return info;
	}

	public:
	using T_ptr = typename View<T, D>::T_ptr;
	/**
	 * @brief Construct a new PatchView with size 0
	 */
	PatchView() : View<T, D + 1>() {}

	/**
	 * @brief Construct a new View object
	 *
	 * @param data pointer to the first element in the patch (non-ghost cell element)
	 * @param strides the strides in each direction
	 * @param lengths the lengths in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @param ldm the local data manager for the data
	 */
	PatchView(T_ptr                         data,
	          const std::array<int, D + 1> &strides,
	          const std::array<int, D + 1> &ghost_start,
	          const std::array<int, D + 1> &start,
	          const std::array<int, D + 1> &end,
	          const std::array<int, D + 1> &ghost_end)
	: View<T, D + 1>(data, strides, ghost_start, start, end, ghost_end)
	{
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
	PatchView(T_ptr data, const std::array<int, D + 1> &strides, const std::array<int, D + 1> &lengths, int num_ghost_cells)
	: View<T, D + 1>(data,
	                 strides,
	                 DetermineGhostStart<D>(num_ghost_cells),
	                 DetermineStart<D>(),
	                 DetermineEnd<D>(lengths),
	                 DetermineGhostEnd<D>(lengths, num_ghost_cells))
	{
	}

	/**
	 * @brief Get the slice on a given face
	 *
	 * @tparam M the dimension of the face
	 * @param f the face
	 * @param offset offset the offset of the value {0,..,0} is the non ghost cell touching the face. {-1,..,-1} is the first ghost cell touching that
	 * face
	 * @return View<T,M+1> a view to the slice on the face
	 */
	template <int M> View<T, M + 1> getSliceOn(Face<D, M> f, const std::array<int, D - M> &offset) const
	{
		SliceInfo<M> info = getSliceOnPriv<M>(f, offset);

		T_ptr new_data = (&(*this)[info.first_value]);
		return View<T, M + 1>(new_data, info.strides, info.ghost_start, info.start, info.end, info.ghost_end);
	}

	/**
	 * @brief Get the gosts slice on a given face
	 *
	 * @tparam M the dimension of the face
	 * @param f the face
	 * @param offset offset the offset of the value {0,..,0} first ghost cell slice touching that face
	 * face
	 * @return View<typename std::remove_const<T>::type, M + 1>
	 */
	template <int M> View<typename std::remove_const<T>::type, M + 1> getGhostSliceOn(Face<D, M> f, const std::array<size_t, D - M> &offset) const
	{
		using noconst_T     = typename std::remove_const<T>::type;
		using noconst_T_ptr = typename std::add_pointer<noconst_T>::type;
		std::array<int, M + 1> new_strides;
		std::array<int, M + 1> new_ghost_start;
		std::array<int, M + 1> new_start;
		std::array<int, M + 1> new_end;
		std::array<int, M + 1> new_ghost_end;
		std::array<int, D + 1> first_value;

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
		new_strides[M]     = this->getStrides()[D];
		new_ghost_start[M] = this->getGhostStart()[D];
		new_start[M]       = this->getStart()[D];
		new_end[M]         = this->getEnd()[D];
		new_ghost_end[M]   = this->getGhostEnd()[D];

		noconst_T_ptr new_data = const_cast<noconst_T_ptr>(&(*this)[first_value]); // Thunderegg doesn't care if values in ghosts are modified
		return View<noconst_T, M + 1>(new_data, new_strides, new_ghost_start, new_start, new_end, new_ghost_end);
	}

	ComponentView<T, D> getComponentView(int component_index) const
	{
		if constexpr (ENABLE_DEBUG) {
			if (component_index < this->getStart()[D] || component_index > this->getEnd()[D]) {
				throw RuntimeError("invalid component index");
			}
		}
		std::array<int, D>     new_strides;
		std::array<int, D>     new_ghost_start;
		std::array<int, D>     new_start;
		std::array<int, D>     new_end;
		std::array<int, D>     new_ghost_end;
		std::array<int, D + 1> first_value;

		first_value    = this->getGhostStart();
		first_value[D] = component_index;

		for (size_t axis = 0; axis < (size_t) D; axis++) {
			new_strides[axis]     = this->getStrides()[axis];
			new_ghost_start[axis] = this->getGhostStart()[axis];
			new_start[axis]       = this->getStart()[axis];
			new_end[axis]         = this->getEnd()[axis];
			new_ghost_end[axis]   = this->getGhostEnd()[axis];
		}

		T_ptr new_data = (&(*this)[first_value]);
		return ComponentView<T, D>(new_data, new_strides, new_ghost_start, new_start, new_end, new_ghost_end);
	}
	operator PatchView<std::add_const_t<T>, D>() const
	{
		return PatchView<std::add_const_t<T>, D>(this->getData() + this->getIndex(this->getGhostStart()),
		                                         this->getStrides(),
		                                         this->getGhostStart(),
		                                         this->getStart(),
		                                         this->getEnd(),
		                                         this->getGhostEnd());
	}
};
extern template class PatchView<double, 1>;
extern template View<double, 1> PatchView<double, 1>::getSliceOn(Face<1, 0>, const std::array<int, 1> &) const;
extern template View<double, 1> PatchView<double, 1>::getGhostSliceOn(Face<1, 0>, const std::array<size_t, 1> &) const;

extern template class PatchView<double, 2>;
extern template View<double, 1> PatchView<double, 2>::getSliceOn(Face<2, 0>, const std::array<int, 2> &) const;
extern template View<double, 2> PatchView<double, 2>::getSliceOn(Face<2, 1>, const std::array<int, 1> &) const;
extern template View<double, 1> PatchView<double, 2>::getGhostSliceOn(Face<2, 0>, const std::array<size_t, 2> &) const;
extern template View<double, 2> PatchView<double, 2>::getGhostSliceOn(Face<2, 1>, const std::array<size_t, 1> &) const;

extern template class PatchView<double, 3>;
extern template View<double, 1> PatchView<double, 3>::getSliceOn(Face<3, 0>, const std::array<int, 3> &) const;
extern template View<double, 2> PatchView<double, 3>::getSliceOn(Face<3, 1>, const std::array<int, 2> &) const;
extern template View<double, 3> PatchView<double, 3>::getSliceOn(Face<3, 2>, const std::array<int, 1> &) const;
extern template View<double, 1> PatchView<double, 3>::getGhostSliceOn(Face<3, 0>, const std::array<size_t, 3> &) const;
extern template View<double, 2> PatchView<double, 3>::getGhostSliceOn(Face<3, 1>, const std::array<size_t, 2> &) const;
extern template View<double, 3> PatchView<double, 3>::getGhostSliceOn(Face<3, 2>, const std::array<size_t, 1> &) const;

extern template class PatchView<const double, 1>;
extern template View<const double, 1> PatchView<const double, 1>::getSliceOn(Face<1, 0>, const std::array<int, 1> &) const;
extern template View<double, 1>       PatchView<const double, 1>::getGhostSliceOn(Face<1, 0>, const std::array<size_t, 1> &) const;

extern template class PatchView<const double, 2>;
extern template View<const double, 1> PatchView<const double, 2>::getSliceOn(Face<2, 0>, const std::array<int, 2> &) const;
extern template View<const double, 2> PatchView<const double, 2>::getSliceOn(Face<2, 1>, const std::array<int, 1> &) const;
extern template View<double, 1>       PatchView<const double, 2>::getGhostSliceOn(Face<2, 0>, const std::array<size_t, 2> &) const;
extern template View<double, 2>       PatchView<const double, 2>::getGhostSliceOn(Face<2, 1>, const std::array<size_t, 1> &) const;

extern template class PatchView<const double, 3>;
extern template View<const double, 1> PatchView<const double, 3>::getSliceOn(Face<3, 0>, const std::array<int, 3> &) const;
extern template View<const double, 2> PatchView<const double, 3>::getSliceOn(Face<3, 1>, const std::array<int, 2> &) const;
extern template View<const double, 3> PatchView<const double, 3>::getSliceOn(Face<3, 2>, const std::array<int, 1> &) const;
extern template View<double, 1>       PatchView<const double, 3>::getGhostSliceOn(Face<3, 0>, const std::array<size_t, 3> &) const;
extern template View<double, 2>       PatchView<const double, 3>::getGhostSliceOn(Face<3, 1>, const std::array<size_t, 2> &) const;
extern template View<double, 3>       PatchView<const double, 3>::getGhostSliceOn(Face<3, 2>, const std::array<size_t, 1> &) const;
} // namespace ThunderEgg
#endif