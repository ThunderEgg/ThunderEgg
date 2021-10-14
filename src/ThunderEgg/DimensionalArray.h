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

#ifndef THUNDEREGG_DIMENSIONALARRAY_H
#define THUNDEREGG_DIMENSIONALARRAY_H
/**
 * @file
 *
 * @brief DimensionalArray class
 */

namespace ThunderEgg
{
/**
 * @brief A dimensional array, used for storing templated with and integer. For example a dimensional array of length 3 holding type T will have three
 * values T<0> at index 0, T<1> at index 1, and T<2> at index 2.
 *
 * @tparam N  the length of the array
 * @tparam T  the type to hold
 */
template <int N, template <int> class T> class DimensionalArray
{
	private:
	/**
	 * @brief The values for the previous indexes
	 *
	 */
	DimensionalArray<N - 1, T> prev;
	/**
	 * @brief The value for index N-1
	 */
	T<N - 1> t;

	public:
	/**
	 * @brief Get a value at a given index
	 *
	 * @tparam I the index
	 * @return T<I>&  the value at that index
	 */
	template <int I> T<I> &get()
	{
		static_assert(I < N, "invalid index value");
		if constexpr (I == N - 1) {
			return t;
		} else {
			return prev.template get<I>();
		}
	}
	/**
	 * @brief Get a value at an index
	 *
	 * @tparam I the index
	 * @return const T<I>& the value
	 */
	template <int I> const T<I> &get() const
	{
		static_assert(I < N, "invalid index value");
		if constexpr (I == N - 1) {
			return t;
		} else {
			return prev.template get<I>();
		}
	}
};
/**
 * @brief Dimensional Array specialization for length of 1
 *
 * @tparam T the type that is being held
 */
template <template <int> class T> class DimensionalArray<1, T>
{
	private:
	/**
	 * @brief The value at index 0
	 */
	T<0> t;

	public:
	/**
	 * @brief Get a value at an index
	 *
	 * @tparam I the index
	 * @return const T<I>& the value
	 */
	template <int I> T<I> &get()
	{
		static_assert(I == 0, "invalid index value");
		return t;
	}
	/**
	 * @brief Get a value at an index
	 *
	 * @tparam I the index
	 * @return const T<I>& the value
	 */
	template <int I> const T<I> &get() const
	{
		static_assert(I == 0, "invalid index value");
		return t;
	}
};
} // namespace ThunderEgg
#endif
