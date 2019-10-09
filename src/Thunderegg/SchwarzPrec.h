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

#ifndef SCHWARZPREC_H
#define SCHWARZPREC_H
#include <Thunderegg/Operator.h>
#include <Thunderegg/Schur/IfaceInterp.h>
#include <Thunderegg/Schur/PatchSolver.h>
#include <Thunderegg/Schur/SchurHelper.h>
namespace Thunderegg
{
/**
 * @brief Additive Schwarz Preconditioner
 */
template <size_t D> class SchwarzPrec : public Operator<D>
{
	private:
	/**
	 * @brief the SchurHelper
	 */
	std::shared_ptr<Schur::SchurHelper<D>> sh;
	std::shared_ptr<Schur::PatchSolver<D>> solver;
	std::shared_ptr<Schur::IfaceInterp<D>> interp;

	public:
	/**
	 * @brief Create new preconditioner
	 *
	 * @param sh the SchurHelper
	 */
	SchwarzPrec(std::shared_ptr<Schur::SchurHelper<D>> sh,
	            std::shared_ptr<Schur::PatchSolver<D>> solver,
	            std::shared_ptr<Schur::IfaceInterp<D>> interp)
	{
		this->sh     = sh;
		this->solver = solver;
		this->interp = interp;
	}
	/**
	 * @brief Apply schwarz preconditioner
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
	{
		auto gamma = sh->getNewSchurDistVec();
		interp->interpolateToInterface(x, gamma);
		solver->solve(x, b, gamma);
	}
};
} // namespace Thunderegg
#endif
