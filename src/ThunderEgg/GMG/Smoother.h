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

#ifndef THUNDEREGG_GMG_SMOOTHER_H
#define THUNDEREGG_GMG_SMOOTHER_H
/**
 * @file
 *
 * @brief Smoother class
 */
#include <ThunderEgg/Vector.h>
namespace ThunderEgg::GMG {
/**
 * @brief Abstract class for smoothing operators.
 */
template<int D>
class Smoother
{
public:
  /**
   * @brief Destroy the Smoother object
   */
  virtual ~Smoother() {}
  /**
   * @brief Clone this smoother
   *
   * @return Smoother<D>* a newly allocated copy of this smoother
   */
  virtual Smoother<D>* clone() const = 0;
  /**
   * @brief Virtual function that derived classes have to implement.
   *
   * @param f the RHS vector
   * @param u the solution vector, updated upon return.
   */
  virtual void smooth(const Vector<D>& f, Vector<D>& u) const = 0;
};
} // namespace ThunderEgg::GMG
#endif