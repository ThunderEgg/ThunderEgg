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

#ifndef THUNDEREGG_FACE_H
#define THUNDEREGG_FACE_H
#include <cstddef>
namespace ThunderEgg
{
template <int D, int M> class Face
{
	private:
	/**
	 * @brief The value of the face
	 */
	unsigned char value;
	/**
	 * @brief 2^N
	 */
	static constexpr size_t twopow(size_t N)
	{
		return 0b1 << N;
	}
	/**
	 * @brief N!
	 */
	static constexpr size_t factorial(size_t N)
	{
		return N == 0 ? 1 : N * factorial(N - 1);
	}
	/**
	 * @brief N choose K
	 */
	static constexpr size_t choose(size_t N, size_t K)
	{
		return factorial(N) / (factorial(K) * factorial(N - K));
	}

	public:
	/**
	 * @brief the number of faces
	 */
	static constexpr size_t num_faces = twopow(D - M) * choose(D, M);
	/**
	 * @brief Construct a new Face object
	 *
	 * @param value the value of the face to give
	 */
	Face(unsigned char value) : value(value)
	{
		static_assert(D <= 3 && D > 0, "Only up to 3 dimensions supported");
		static_assert(M > 0 && M < D, "Invalid M value");
	}
	size_t getIndex()
	{
		return value;
	}
};

} // namespace ThunderEgg
#endif
