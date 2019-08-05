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

#ifndef THUNDEREGG_DOMAIN_TOOLS_H
#define THUNDEREGG_DOMAIN_TOOLS_H

#include <Thunderegg/Domain.h>
#include <functional>

namespace Thunderegg
{
template <size_t D> struct DomainTools {
	static void getRealCoord(std::shared_ptr<PatchInfo<D>> pinfo, const std::array<int, D> &coord,
	                         std::array<double, D> &real_coord)
	{
		loop<0, D - 1>([&](int dir) {
			if (coord[dir] == -1) {
				real_coord[dir] = pinfo->starts[dir];
			} else if (coord[dir] == pinfo->ns[dir]) {
				real_coord[dir] = pinfo->starts[dir] + pinfo->spacings[dir] * pinfo->ns[dir];
			} else {
				real_coord[dir] = pinfo->starts[dir] + pinfo->spacings[dir] / 2.0
				                  + pinfo->spacings[dir] * coord[dir];
			}
		});
	}
	static void getRealCoordBound(std::shared_ptr<PatchInfo<D>> pinfo,
	                              const std::array<int, D - 1> &coord, Side<D> s,
	                              std::array<double, D> &real_coord)
	{
		for (size_t dir = 0; dir < s.axis(); dir++) {
			if (coord[dir] == -1) {
				real_coord[dir] = pinfo->starts[dir];
			} else if (coord[dir] == pinfo->ns[dir]) {
				real_coord[dir] = pinfo->starts[dir] + pinfo->spacings[dir] * pinfo->ns[dir];
			} else {
				real_coord[dir] = pinfo->starts[dir] + pinfo->spacings[dir] / 2.0
				                  + pinfo->spacings[dir] * coord[dir];
			}
		}
		if (s.isLowerOnAxis()) {
			real_coord[s.axis()] = pinfo->starts[s.axis()];
		} else {
			real_coord[s.axis()]
			= pinfo->starts[s.axis()] + pinfo->spacings[s.axis()] * pinfo->ns[s.axis()];
		}
		for (size_t dir = s.axis() + 1; dir < D; dir++) {
			if (coord[dir - 1] == -1) {
				real_coord[dir] = pinfo->starts[dir];
			} else if (coord[dir - 1] == pinfo->ns[dir]) {
				real_coord[dir] = pinfo->starts[dir] + pinfo->spacings[dir] * pinfo->ns[dir];
			} else {
				real_coord[dir] = pinfo->starts[dir] + pinfo->spacings[dir] / 2.0
				                  + pinfo->spacings[dir] * coord[dir - 1];
			}
		}
	}
	static void setValues(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec,
	                      std::function<double(const std::array<double, D> &)> func)
	{
		std::array<double, D> real_coord;
		for (int i = 0; i < vec->getNumLocalPatches(); i++) {
			LocalData<D> ld    = vec->getLocalData(i);
			auto         pinfo = domain->getPatchInfoVector()[i];
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
				getRealCoord(pinfo, coord, real_coord);
				ld[coord] = func(real_coord);
			});
		}
	}
	static void setBCValues(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D - 1>> vec,
	                        std::function<double(const std::array<double, D> &)> func)
	{
		std::array<double, D> real_coord;
		for (int i = 0; i < vec->getNumLocalPatches(); i++) {
			auto pinfo = domain->getPatchInfoMap()[domain->patch_id_bc_map_vec[i]];
			for (Side<D> s : Side<D>::getValues()) {
				if (!pinfo->hasNbr(s)) {
					LocalData<D - 1> ld = vec->getLocalData(pinfo->getBCLocalIndex(s));
					nested_loop<D - 1>(ld.getStart(), ld.getEnd(),
					                   [&](const std::array<int, D - 1> &coord) {
						                   getRealCoordBound(pinfo, coord, s, real_coord);
						                   ld[coord] = func(real_coord);
					                   });
				}
			}
		}
	}
};
} // namespace Thunderegg
#endif