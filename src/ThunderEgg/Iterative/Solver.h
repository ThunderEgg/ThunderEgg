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

namespace ThunderEgg
{
namespace Iterative
{
template <int D> class Solver
{
	public:
	virtual int solve(std::shared_ptr<VectorGenerator<D>> vg, std::shared_ptr<const Operator<D>> A,
	                  std::shared_ptr<Vector<D>> x, std::shared_ptr<const Vector<D>> b,
	                  std::shared_ptr<const Operator<D>> Mr = nullptr)
	= 0;
};
//
} // namespace Iterative
} // namespace ThunderEgg