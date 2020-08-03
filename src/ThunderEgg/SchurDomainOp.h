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

#ifndef THUNDEREGG_SCHURDOMAINOP_H
#define THUNDEREGG_SCHURDOMAINOP_H

#include <ThunderEgg/Operator.h>
#include <ThunderEgg/Schur/PatchOperator.h>
#include <ThunderEgg/Schur/SchurHelper.h>

namespace ThunderEgg
{
template <size_t D> class SchurDomainOp : public Operator<D>
{
	private:
	/**
	 * @brief PETSc Matrix object
	 */
	std::shared_ptr<Schur::SchurHelper<D>>   helper;
	std::shared_ptr<Schur::IfaceInterp<D>>   interp;
	std::shared_ptr<Schur::PatchOperator<D>> op;

	public:
	/**
	 * @brief Crate new WrapOp
	 *
	 * @param matrix the PETSc matrix
	 */
	SchurDomainOp(std::shared_ptr<Schur::SchurHelper<D>>   helper,
	              std::shared_ptr<Schur::IfaceInterp<D>>   interp,
	              std::shared_ptr<Schur::PatchOperator<D>> op)
	{
		this->helper = helper;
		this->interp = interp;
		this->op     = op;
	}
	/**
	 * @brief Perform matrix/vector multiply.
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
	{
		auto gamma = helper->getNewSchurDistVec();
		interp->interpolateToInterface(x, gamma);
		op->apply(x, gamma, b);
	}
};
} // namespace ThunderEgg
#endif
