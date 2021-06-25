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

#ifndef THUNDEREGG_ITERATIVE_SOLVER_H
#define THUNDEREGG_ITERATIVE_SOLVER_H

#include <ThunderEgg/Operator.h>
#include <ThunderEgg/VectorGenerator.h>
#include <iostream>

namespace ThunderEgg
{
namespace Iterative
{
/**
 * @brief Represents an iterative solver
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D> class Solver
{
	public:
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
	virtual int solve(const Operator<D> &A,
	                  Vector<D> &        x,
	                  const Vector<D> &  b,
	                  const Operator<D> *Mr     = nullptr,
	                  bool               output = false,
	                  std::ostream &     os     = std::cout) const = 0;
};
//
} // namespace Iterative
} // namespace ThunderEgg
#endif