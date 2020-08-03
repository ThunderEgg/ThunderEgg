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

#ifndef THUNDEREGG_SCHUR_IFACEINTERP_H
#define THUNDEREGG_SCHUR_IFACEINTERP_H
#include <ThunderEgg/GMG/CycleFactoryCtx.h>
#include <ThunderEgg/Schur/IfaceType.h>
#include <ThunderEgg/Vector.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief An abstract class that interpolates to the interfaces in the Schur compliment system.
 *
 * @tparam D the number of Cartesian dimensions in the patches.
 */
template <size_t D> class IfaceInterp
{
	public:
	virtual ~IfaceInterp() {}

	/**
	 * @brief Interpolate the vector to the interfaces
	 *
	 * @param u the input vector
	 * @param interp the interface vector to be interpolated to
	 */
	virtual void interpolateToInterface(std::shared_ptr<const Vector<D>> u,
	                                    std::shared_ptr<Vector<D - 1>>   interp)
	= 0;
	virtual std::shared_ptr<IfaceInterp<D>> getNewIfaceInterp(GMG::CycleFactoryCtx<D> ctx) = 0;
};
} // namespace Schur
} // namespace ThunderEgg
#endif
