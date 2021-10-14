/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2018-2021 Scott Aiton
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

#ifndef THUNDEREGG_GMG_INTERPOLATOR_H
#define THUNDEREGG_GMG_INTERPOLATOR_H
/**
 * @file
 *
 * @brief Interpolator class
 */
#include <ThunderEgg/Vector.h>

namespace ThunderEgg::GMG
{
/**
 * @brief Abstract class for interpolation operators.
 */
template <int D> class Interpolator
{
	public:
	/**
	 * @brief Destroy the Interpolator object
	 */
	virtual ~Interpolator() {}
	/**
	 * @brief Clone this interpolator
	 *
	 * @return Interpolator<D>* a newly allocated copy of this interpolator
	 */
	virtual Interpolator<D> *clone() const = 0;
	/**
	 * @brief Virtual interpolation operation that needs to be implemented in derived classes.
	 *
	 * @param coarse the input vector from the coarser level.
	 * @param fine the output vector for the fine level.
	 */
	virtual void interpolate(const Vector<D> &coarse, Vector<D> &fine) const = 0;
};
} // namespace ThunderEgg::GMG
#endif