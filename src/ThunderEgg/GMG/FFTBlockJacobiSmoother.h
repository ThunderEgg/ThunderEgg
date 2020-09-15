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

#ifndef THUNDEREGG_GMG_BLOCKJACOBISMOOTHER_H
#define THUNDEREGG_GMG_BLOCKJACOBISMOOTHER_H

#include <ThunderEgg/Schur/InterfaceDomain.h>

namespace ThunderEgg
{
namespace GMG
{
/**
 * @brief A block Jacobi smoother that uses FFTW solves on each patch. Implemented using the
 * InterfaceDomain class.
 */
template <int D> class FFTBlockJacobiSmoother : public Smoother<D>
{
	private:
	/**
	 * @brief point to the InterfaceDomain object.
	 */
	std::shared_ptr<Schur::InterfaceDomain<D>> sh;
	std::shared_ptr<Schur::PatchSolver<D>>     solver;
	std::shared_ptr<Schur::IfaceInterp<D>>     interp;

	public:
	/**
	 * @brief Create new smoother with InterfaceDomain object
	 *
	 * @param sh pointer to the InterfaceDomain object
	 */
	FFTBlockJacobiSmoother(std::shared_ptr<Schur::InterfaceDomain<D>> sh,
	                       std::shared_ptr<Schur::PatchSolver<D>>     solver,
	                       std::shared_ptr<Schur::IfaceInterp<D>>     interp)
	{
		this->sh     = sh;
		this->solver = solver;
		this->interp = interp;
	}
	/**
	 * @brief Run an iteration of smoothing.
	 *
	 * @param f the RHS vector
	 * @param u the solution vector, updated upon return.
	 */
	void smooth(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u) const
	{
		auto gamma = sh->getNewPatchLocalInterfaceVector();
		interp->interpolateToInterface(u, gamma);
		solver->solve(f, u, gamma);
	}
};
} // namespace GMG
} // namespace ThunderEgg
#endif