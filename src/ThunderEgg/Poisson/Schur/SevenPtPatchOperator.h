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

#ifndef THUNDEREGG_POISSON_SCHUR_SEVENPTPATCHOPERATOR_H
#define THUNDEREGG_POISSON_SCHUR_SEVENPTPATCHOPERATOR_H

#include <ThunderEgg/Schur/PatchOperator.h>

namespace ThunderEgg
{
namespace Poisson
{
namespace Schur
{
class SevenPtPatchOperator : public ThunderEgg::Schur::PatchOperator<3>
{
	public:
	void apply(std::shared_ptr<const Vector<3>> u, std::shared_ptr<const Vector<2>> gamma,
	           std::shared_ptr<Vector<3>> f) override;
	void applyWithInterface(ThunderEgg::Schur::SchurInfo<3> &d, const LocalData<3> u,
	                        std::shared_ptr<const Vector<2>> gamma, LocalData<3> f) override;
	void addInterfaceToRHS(ThunderEgg::Schur::SchurInfo<3> &sinfo,
	                       std::shared_ptr<const Vector<2>> gamma, LocalData<3> f) override;
	void apply(const ThunderEgg::Schur::SchurInfo<3> &sinfo, const LocalData<3> u,
	           LocalData<3> f) override;
	std::shared_ptr<ThunderEgg::Schur::PatchOperator<3>>
	getNewPatchOperator(GMG::CycleFactoryCtx<3> ctx) override
	{
		return std::shared_ptr<PatchOperator<3>>(new SevenPtPatchOperator());
	}
};
} // namespace Schur
} // namespace Poisson
} // namespace ThunderEgg
#endif