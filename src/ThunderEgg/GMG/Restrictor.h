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

#ifndef THUNDEREGG_GMG_RESTRICTOR_H
#define THUNDEREGG_GMG_RESTRICTOR_H

#include <ThunderEgg/Vector.h>

namespace ThunderEgg
{
namespace GMG
{
/**
 * @brief Base class for multi-grid restriction operators.
 */
template <int D> class Restrictor
{
	public:
	/**
	 * @brief Destroy the Restrictor object
	 */
	virtual ~Restrictor() {}
	/**
	 * @brief Virtual function that base classes have to implement.
	 *
	 * @param fine the input vector that is restricted.
	 * @param coarse the output vector that is restricted to.
	 */

	/**
	 * @brief
	 *
	 * @param fine
	 */
	virtual Vector<D> restrict(const Vector<D> &fine) const = 0;
	/**
	 * @brief Get get a new vector for the coarser domain
	 *
	 * @param num_components the number of components in the vector
	 * @return Vector<D> the vector
	 */
	virtual Vector<D> getNewCoarserVector(int num_components) const = 0;
};
} // namespace GMG
} // namespace ThunderEgg
#endif