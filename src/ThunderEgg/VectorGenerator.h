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

#ifndef THUNDEREGG_VECTORGENERATOR_H
#define THUNDEREGG_VECTORGENERATOR_H
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/Side.h>
#include <ThunderEgg/Vector.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <numeric>
namespace ThunderEgg
{
/**
 * @brief Generates temporary work vectors
 *
 * This is used in various ThunderEgg classes to generate needed work vectors
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class VectorGenerator
{
	public:
	/**
	 * @brief Destroy the VectorGenerator object
	 */
	virtual ~VectorGenerator() {}
	/**
	 * @brief Get a new Vector
	 *
	 * @return std::shared_ptr<Vector<D>> the Vector
	 */
	virtual std::shared_ptr<Vector<D>> getNewVector() = 0;
};
} // namespace ThunderEgg
#endif
