/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_LOOPS_H
#define THUNDEREGG_LOOPS_H

#include <array>

namespace Thunderegg
{
template <int start, int stop, typename T> class Loop
{
	public:
	static void inline loop_loop(T lambda)
	{
		lambda(start);
		Loop<start + 1, stop, T>::loop_loop(lambda);
	}
};

template <int start, typename T> class Loop<start, start, T>
{
	public:
	static void inline loop_loop(T lambda)
	{
		lambda(start);
	}
};
template <int start, int stop, typename T> inline void loop(T lambda)
{
	Loop<start, stop, T>::loop_loop(lambda);
}
template <size_t D, size_t Dir, typename T> class NestedLoop
{
	public:
	static void inline nested_loop_loop(std::array<int, D> &coord, std::array<int, D> &start,
	                                    std::array<int, D> &end, T lambda)
	{
		for (coord[Dir] = start[Dir]; coord[Dir] <= end[Dir]; coord[Dir]++) {
			NestedLoop<D, Dir - 1, T>::nested_loop_loop(coord, start, end, lambda);
		}
	}
};

template <size_t D, typename T> class NestedLoop<D, 0, T>
{
	public:
	static void inline nested_loop_loop(std::array<int, D> &coord, std::array<int, D> &start,
	                                    std::array<int, D> &end, T lambda)
	{
		for (coord[0] = start[0]; coord[0] <= end[0]; coord[0]++) {
			lambda(coord);
		}
	}
};
template <size_t D, typename T>
inline void nested_loop(std::array<int, D> start, std::array<int, D> end, T lambda)
{
	std::array<int, D> coord = start;
	NestedLoop<D, D - 1, T>::nested_loop_loop(coord, start, end, lambda);
}
} // namespace Thunderegg
#endif