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

#ifndef THUNDEREGG_SCHUR_PATCHSOLVERWRAPPER_H
#define THUNDEREGG_SCHUR_PATCHSOLVERWRAPPER_H

#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>

namespace ThunderEgg
{
namespace Schur
{
template <int D> class PatchSolverWrapper : public Operator<D - 1>
{
	private:
	std::shared_ptr<InterfaceDomain<D>> iface_domain;
	std::shared_ptr<PatchSolver<D>>     solver;

	public:
	PatchSolverWrapper(std::shared_ptr<InterfaceDomain<D>> iface_domain,
	                   std::shared_ptr<PatchSolver<D>>     solver)
	{
		this->iface_domain = iface_domain;
		this->solver       = solver;
	}
	/**
	 * @brief Apply Schur matrix
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D - 1>> x,
	           std::shared_ptr<Vector<D - 1>>       b) const override
	{
	}
	/**
	 * @brief Get the RHS for the Schur system from a given RHS for the domain system
	 *
	 * @param domain_b the domain rhs
	 * @param schur_b the Schur rhs
	 */
	void getSchurRHSFromDomainRHS(std::shared_ptr<const Vector<D>> domain_b,
	                              std::shared_ptr<Vector<D - 1>>   schur_b) const
	{
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
