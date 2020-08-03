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

#ifndef THUNDEREGG_SCHUR_BILINEARINTERPOLATOR_H
#define THUNDEREGG_SCHUR_BILINEARINTERPOLATOR_H

#include <ThunderEgg/Schur/IfaceInterp.h>
#include <ThunderEgg/Schur/SchurHelper.h>

namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief a bilinear interpolator for interface values
 */
class BilinearInterpolator : public IfaceInterp<2>
{
	private:
	std::shared_ptr<SchurHelper<2>> sh;

	public:
	BilinearInterpolator(std::shared_ptr<SchurHelper<2>> sh)
	{
		this->sh = sh;
	}
	void interpolateToInterface(std::shared_ptr<const Vector<2>> u,
	                            std::shared_ptr<Vector<1>>       interp) override;
	void interpolate(const std::vector<std::shared_ptr<SchurInfo<2>>> &patches,
	                 std::shared_ptr<const Vector<2>> u, std::shared_ptr<Vector<1>> interp);
	void interpolate(SchurInfo<2> &d, std::shared_ptr<const Vector<2>> u,
	                 std::shared_ptr<Vector<1>> interp);
	void interpolate(SchurInfo<2> &d, Side<2> s, int local_index, IfaceType<2> itype,
	                 std::shared_ptr<const Vector<2>> u, std::shared_ptr<Vector<1>> interp);
	std::shared_ptr<IfaceInterp<2>> getNewIfaceInterp(GMG::CycleFactoryCtx<2> ctx) override
	{
		return std::shared_ptr<IfaceInterp<2>>(new BilinearInterpolator(ctx.sh));
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
