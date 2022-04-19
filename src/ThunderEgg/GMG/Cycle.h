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

#ifndef THUNDEREGG_GMG_CYCLE_H
#define THUNDEREGG_GMG_CYCLE_H
/**
 * @file
 *
 * @brief Cycle class
 */

#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/Vector.h>
#include <list>

namespace ThunderEgg::GMG {

/**
 * @brief Base abstract class for cycles.
 *
 * There is an abstract visit() function for base classes to implement.
 */
template<int D>
  requires is_supported_dimension<D>
class Cycle : public Operator<D>
{
private:
  /**
   * @brief Implimentation class
   */
  class Implimentation;

  /**
   * @brief pointer to the implimentation
   */
  std::shared_ptr<const Implimentation> implimentation;

protected:
  /**
   * @brief Prepare vectors for coarser level.
   *
   * @param level the current level
   * @param f the rhs vector cooresponding to the level
   * @param u the solution vector cooresponding to the level
   * @return Vector<D> the restricted residual vector
   */
  Vector<D> restrict(const Level<D>& level, const Vector<D>& f, const Vector<D>& u) const;

  /**
   * @brief Virtual visit function that needs to be implemented in derived classes.
   *
   * @param level the level currently begin visited.
   * @param f the rhs vector cooresponding to the level
   * @param u the solution vector cooresponding to the level
   */
  virtual void
  visit(const Level<D>& level, const Vector<D>& f, Vector<D>& u) const = 0;

public:
  /**
   * @brief Create new cycle object.
   *
   * @param finest_level the finest level object.
   */
  Cycle(const Level<D>& finest_level);

  /**
   * @brief Run one iteration of the cycle.
   *
   * Performs one cycle on on the system `Au=f` where `A` is the operator for the
   * finest level.
   *
   * @param f the RHS vector.
   * @param u the solution vector.
   */
  void
  apply(const Vector<D>& f, Vector<D>& u) const override;

  /**
   * @brief Get the finest Level
   *
   * @return const Level<D>& the Level
   */
  const Level<D>&
  getFinestLevel() const;
};
extern template class Cycle<2>;
extern template class Cycle<3>;
} // namespace ThunderEgg::GMG
#endif