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

#ifndef THUNDEREGG_GMG_CYCLEFACTORYCTX_H
#define THUNDEREGG_GMG_CYCLEFACTORYCTX_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/GMG/Level.h>
#include <Thunderegg/SchurHelper.h>
namespace Thunderegg
{
template<size_t D>
class PatchOperator;
namespace GMG
{
template <size_t D> struct CycleFactoryCtx {
	std::shared_ptr<Domain<D>>      domain;
	std::shared_ptr<SchurHelper<D>> sh;
	std::shared_ptr<Level<D>>       level;
	std::shared_ptr<PatchOperator<D>>       op;
};
} // namespace GMG
} // namespace Thunderegg
#endif