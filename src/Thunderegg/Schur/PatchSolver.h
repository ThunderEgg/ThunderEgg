/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_SCHUR_PATCHSOLVER_H
#define THUNDEREGG_SCHUR_PATCHSOLVER_H

#include <Thunderegg/GMG/CycleFactoryCtx.h>
#include <Thunderegg/Schur/SchurInfo.h>
#include <Thunderegg/Vector.h>

namespace Thunderegg
{
namespace Schur
{
/**
 * @brief Solves the problem on the patches using a specified interface value
 *
 * @tparam D the number of cartesian dimensions
 */
template <size_t D> class PatchSolver
{
	public:
	/**
	 * @brief Destroy the Patch Solver object
	 */
	virtual ~PatchSolver() {}
	/**
	 * @brief add a patch to the solver
	 *
	 * The solver do any necessary setup for the patch
	 *
	 * @param sinfo the patch
	 */
	virtual void addPatch(SchurInfo<D> &sinfo) = 0;
	/**
	 * @brief Solve all the patches in the domain
	 *
	 * @param patches the patches
	 * @param f the rhs vector
	 * @param u the lhs vector
	 * @param gamma the interface values to use
	 */
	virtual void solve(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u,
	                   std::shared_ptr<const Vector<D - 1>> gamma)
	= 0;
	virtual std::shared_ptr<PatchSolver<D>> getNewPatchSolver(GMG::CycleFactoryCtx<D> ctx) = 0;
};
} // namespace Schur
} // namespace Thunderegg
#endif
