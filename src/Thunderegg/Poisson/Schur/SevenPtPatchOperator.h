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

#ifndef THUNDEREGG_POISSON_SCHUR_SEVENPTPATCHOPERATOR_H
#define THUNDEREGG_POISSON_SCHUR_SEVENPTPATCHOPERATOR_H

#include <Thunderegg/Schur/PatchOperator.h>

namespace Thunderegg
{
namespace Poisson
{
namespace Schur
{
class SevenPtPatchOperator : public Thunderegg::Schur::PatchOperator<3>
{
	public:
	void apply(std::shared_ptr<const Vector<3>> u, std::shared_ptr<const Vector<2>> gamma,
	           std::shared_ptr<Vector<3>> f) override;
	void applyWithInterface(Thunderegg::Schur::SchurInfo<3> &d, const LocalData<3> u,
	                        std::shared_ptr<const Vector<2>> gamma, LocalData<3> f);
	void addInterfaceToRHS(Thunderegg::Schur::SchurInfo<3> &sinfo,
	                       std::shared_ptr<const Vector<2>> gamma, LocalData<3> f);
	void apply(const Thunderegg::Schur::SchurInfo<3> &sinfo, const LocalData<3> u, LocalData<3> f);
	std::shared_ptr<Thunderegg::Schur::PatchOperator<3>>
	getNewPatchOperator(GMG::CycleFactoryCtx<3> ctx)
	{
		return std::shared_ptr<PatchOperator<3>>(new SevenPtPatchOperator());
	}
};
} // namespace Schur
} // namespace Poisson
} // namespace Thunderegg
#endif
