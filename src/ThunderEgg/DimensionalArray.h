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

#ifndef THUNDEREGG_DIMENSIONALARRAY_H
#define THUNDEREGG_DIMENSIONALARRAY_H

namespace ThunderEgg
{
template <int N, template <int> class T> class DimensionalArray
{
	private:
	T<N>                       t;
	DimensionalArray<N - 1, T> next;

	public:
	template <int I> T<I> &get()
	{
		static_assert(I < N, "invalid index value");
		return next.template get<I>();
	}
	template <> T<N> &get()
	{
		return t;
	}
	template <int I> const T<I> &get() const
	{
		static_assert(I < N, "invalid index value");
		return next.template get<I>();
	}
	template <> const T<N> &get() const
	{
		return t;
	}
};
template <template <int> class T> class DimensionalArray<0, T>
{
	private:
	T<0> t;

	public:
	template <int I> T<I> &get()
	{
		static_assert(I == 0, "invalid index value");
	}
	template <> T<0> &get()
	{
		return t;
	}
	template <int I> const T<I> &get() const
	{
		static_assert(I == 0, "invalid index value");
	}
	template <> const T<0> &get() const
	{
		return t;
	}
};
} // namespace ThunderEgg
#endif
