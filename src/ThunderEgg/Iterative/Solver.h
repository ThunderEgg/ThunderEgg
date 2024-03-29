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

#ifndef THUNDEREGG_ITERATIVE_SOLVER_H
#define THUNDEREGG_ITERATIVE_SOLVER_H
/**
 * @file
 *
 * @brief Solver class
 */

#include <ThunderEgg/Operator.h>
#include <iostream>

namespace ThunderEgg::Iterative {
/**
 * @brief Abstract interface for Iterative solvers
 *
 * @tparam D the number of cartesian dimensions
 */
template<int D>
class Solver
{
public:
  /**
   * @brief Clone this solver
   *
   * @return Solver* a newly allocated copy of this solver
   */
  virtual Solver<D>* clone() const = 0;
  /**
   * @brief Destroy the Solver object
   */
  virtual ~Solver() {}
  /**
   * @brief Perform an iterative solve
   *
   * @param A the matrix
   * @param x the initial LHS guess.
   * @param b the RHS vector.
   * @param Mr the right preconditioner. Set to nullptr if there is no right preconditioner.
   * @param output print output to the provided stream
   * @param os the stream to output to
   *
   * @return the number of iterations
   */
  virtual int solve(const Operator<D>& A,
                    Vector<D>& x,
                    const Vector<D>& b,
                    const Operator<D>* Mr = nullptr,
                    bool output = false,
                    std::ostream& os = std::cout) const = 0;
};
} // namespace ThunderEgg::Iterative
#endif