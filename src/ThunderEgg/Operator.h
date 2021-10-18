/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2018-2021 Scott Aiton
 *  Copyright (c) 2020      Donna Calhoun
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

#ifndef THUNDEREGG_OPERATOR_H
#define THUNDEREGG_OPERATOR_H
/**
 * @file
 *
 * @brief Operator class
 */

#include <ThunderEgg/Vector.h>

namespace ThunderEgg {
/**
 * @brief Base class for operators
 */
template<int D>
class Operator
{
public:
  /**
   * @brief Destroy the Operator object
   */
  virtual ~Operator() {}

  /**
   * @brief Clone this operator
   *
   * @return Operator<D>* a newly allocated copy of this operator
   */
  virtual Operator<D>* clone() const = 0;

  /**
   * @brief Virtual function that base classes have to implement.
   *
   * @param x the input vector.
   * @param b the output vector.
   */
  virtual void apply(const Vector<D>& x, Vector<D>& b) const = 0;
};
} // namespace ThunderEgg
#endif
