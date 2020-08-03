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

#ifndef THUNDEREGG_DOMAINWRAPOP_H
#define THUNDEREGG_DOMAINWRAPOP_H

#include <ThunderEgg/Operator.h>
#include <ThunderEgg/Schur/IfaceInterp.h>
#include <ThunderEgg/Schur/PatchOperator.h>
#include <ThunderEgg/Schur/SchurHelper.h>

namespace ThunderEgg
{
/**
 * @brief Base class for operators
 */
template <size_t D> class DomainWrapOp : public Operator<D>
{
	private:
	std::shared_ptr<Schur::SchurHelper<D>>   sh;
	std::shared_ptr<Schur::PatchOperator<D>> op;
	std::shared_ptr<Schur::IfaceInterp<D>>   interp;

	public:
	DomainWrapOp(std::shared_ptr<Schur::SchurHelper<D>>   sh,
	             std::shared_ptr<Schur::IfaceInterp<D>>   interp,
	             std::shared_ptr<Schur::PatchOperator<D>> op)
	{
		this->sh     = sh;
		this->interp = interp;
		this->op     = op;
	}
	/**
	 * @brief Apply Schur matrix
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
	{
		auto gamma = sh->getNewSchurDistVec();
		interp->interpolateToInterface(x, gamma);
		op->apply(x, gamma, b);
	}
};
} // namespace ThunderEgg
#endif
