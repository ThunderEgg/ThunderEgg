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

#ifndef THUNDEREGG_GMG_CYCLEFACTORYCTX_H
#define THUNDEREGG_GMG_CYCLEFACTORYCTX_H
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/Schur/SchurHelper.h>
namespace ThunderEgg
{
namespace Schur
{
template <size_t D> class PatchOperator;
template <size_t D> class IfaceInterp;
} // namespace Schur
namespace GMG
{
template <size_t D> struct CycleFactoryCtx {
	std::shared_ptr<Domain<D>>               domain;
	std::shared_ptr<Schur::SchurHelper<D>>   sh;
	std::shared_ptr<Level<D>>                coarser_level;
	std::shared_ptr<Level<D>>                finer_level;
	std::shared_ptr<Schur::PatchOperator<D>> op;
	std::shared_ptr<Schur::IfaceInterp<D>>   interp;
};
} // namespace GMG
} // namespace ThunderEgg
#endif