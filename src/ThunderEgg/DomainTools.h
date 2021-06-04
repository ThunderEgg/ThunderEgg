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

#ifndef THUNDEREGG_DOMAIN_TOOLS_H
#define THUNDEREGG_DOMAIN_TOOLS_H

#include <ThunderEgg/Domain.h>
#include <ThunderEgg/RuntimeError.h>

namespace ThunderEgg
{
/**
 * @brief Various tools for filling in values in a domain
 */
namespace DomainTools
{
/**
 * @brief Given a path info object, get the coordinate from a given index into the patch.
 *
 * If one of the values is -1 or ns[axis] it will give the coordinate of the cooresponding
 * interface
 *
 * @tparam D the number of cartesian dimensions
 * @param pinfo the patch to get the coordinate for
 * @param coord the index in the patch
 * @param real_coord (output) the coordnitate of the index
 */
template <int D> void GetRealCoord(const PatchInfo<D> &pinfo, const std::array<int, D> &coord, std::array<double, D> &real_coord)
{
	loop<0, D - 1>([&](int dir) {
		if (coord[dir] == -1) {
			real_coord[dir] = pinfo.starts[dir];
		} else if (coord[dir] == pinfo.ns[dir]) {
			real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] * pinfo.ns[dir];
		} else {
			real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir];
		}
	});
}
/**
 * @brief Given a path info object, get the coordinate from a given index into the patch.
 *
 * @tparam D the number of cartesian dimensions
 * @param pinfo the patch to get the coordinate for
 * @param coord the index in the patch
 * @param real_coord (output) the coordnitate of the index
 */
template <int D> void GetRealCoordGhost(const PatchInfo<D> &pinfo, const std::array<int, D> &coord, std::array<double, D> &real_coord)
{
	loop<0, D - 1>([&](int dir) { real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir]; });
}
/**
 * @brief Given a path info object and a side of the patch, get the coordinate from a given
 * index into the interface of the patch.
 *
 * @tparam D the number of cartesian dimensions
 * @param pinfo the patch to get the coordinate for
 * @param coord the index in the patch
 * @param s the side of the patch that the boundary is on
 * @param real_coord (output) the coordnitate of the index
 */
template <int D> void GetRealCoordBound(const PatchInfo<D> &pinfo, const std::array<int, D - 1> &coord, Side<D> s, std::array<double, D> &real_coord)
{
	for (size_t dir = 0; dir < s.getAxisIndex(); dir++) {
		if (coord[dir] == -1) {
			real_coord[dir] = pinfo.starts[dir];
		} else if (coord[dir] == pinfo.ns[dir]) {
			real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] * pinfo.ns[dir];
		} else {
			real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir];
		}
	}
	if (s.isLowerOnAxis()) {
		real_coord[s.getAxisIndex()] = pinfo.starts[s.getAxisIndex()];
	} else {
		real_coord[s.getAxisIndex()] = pinfo.starts[s.getAxisIndex()] + pinfo.spacings[s.getAxisIndex()] * pinfo.ns[s.getAxisIndex()];
	}
	for (size_t dir = s.getAxisIndex() + 1; dir < D; dir++) {
		if (coord[dir - 1] == -1) {
			real_coord[dir] = pinfo.starts[dir];
		} else if (coord[dir - 1] == pinfo.ns[dir]) {
			real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] * pinfo.ns[dir];
		} else {
			real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir - 1];
		}
	}
}
/**
 * @brief Set the values for a vector with the given function.
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 */
template <int D, typename T> void SetValues(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, int component_index, T func)
{
	if (component_index >= vec->getNumComponents()) {
		throw RuntimeError("Invalid component to set");
	}
	std::array<double, D> real_coord;
	for (int i = 0; i < vec->getNumLocalPatches(); i++) {
		ComponentView<D> ld    = vec->getComponentView(component_index, i);
		auto             pinfo = domain->getPatchInfoVector()[i];
		nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
			GetRealCoord<D>(pinfo, coord, real_coord);
			ld[coord] = func(real_coord);
		});
	}
}
/**
 * @brief This is not intended to be called directly. Set the values for a vector with the given
 * function.
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 */
template <int D, typename T> void _SetValues(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, int component_index, T func)
{
	SetValues(domain, vec, component_index, func);
}
/**
 * @brief This is not intended to be called directly. Set the values for a vector with the given
 * function.
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @tparam Args additional functions
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 * @param args additional functions for additional components
 */
template <int D, typename T, typename... Args>
void _SetValues(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, int component_index, T func, Args... args)
{
	SetValues(domain, vec, component_index, func);
	_SetValues(domain, vec, component_index + 1, args...);
}
/**
 * @brief Set the values for a vector with the given functions
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @tparam Args additional functions
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 * @param args additional functions for additional components
 */
template <int D, typename T, typename... Args> void SetValues(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, T func, Args... args)
{
	_SetValues(domain, vec, 0, func, args...);
}
/**
 * @brief Set the values (including ghost values) for a vector with the given function.
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 */
template <int D, typename T> void SetValuesWithGhost(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, int component_index, T func)
{
	if (component_index >= vec->getNumComponents()) {
		throw RuntimeError("Invalid component to set");
	}
	std::array<double, D> real_coord;
	for (int i = 0; i < vec->getNumLocalPatches(); i++) {
		ComponentView<D> ld    = vec->getComponentView(component_index, i);
		auto             pinfo = domain->getPatchInfoVector()[i];
		nested_loop<D>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, D> &coord) {
			GetRealCoordGhost<D>(pinfo, coord, real_coord);
			ld[coord] = func(real_coord);
		});
	}
}
/**
 * @brief This is not intended to be called directly. Set the values (including ghost values) for a
 * vector with the given function.
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 */
template <int D, typename T> void _SetValuesWithGhost(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, int component_index, T func)
{
	SetValuesWithGhost(domain, vec, component_index, func);
}
/**
 * @brief This is not intended to be called directly. Set the values (including ghost values) for a
 * vector with the given function.
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @tparam Args additional functions
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 * @param args additional functions for additional components
 */
template <int D, typename T, typename... Args>
void _SetValuesWithGhost(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, int component_index, T func, Args... args)
{
	SetValuesWithGhost(domain, vec, component_index, func);
	_SetValuesWithGhost(domain, vec, component_index + 1, args...);
}
/**
 * @brief Set the values (including ghost values) for a vector with the given functions
 *
 * @tparam D the number of cartesian dimensions
 * @tparam T function type
 * @tparam Args additional functions
 * @param domain the domain that we are setting values for
 * @param vec the vector to set values in
 * @param component_index the component to set
 * @param func the function
 * @param args additional functions for additional components
 */
template <int D, typename T, typename... Args>
void SetValuesWithGhost(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> vec, T func, Args... args)
{
	_SetValuesWithGhost(domain, vec, 0, func, args...);
}
/**
 * @brief Set the value of a boundary vector using a given function.
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D, typename T> void SetBCValues(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D - 1>> vec, T func, int component_index = 0)
{
	if (component_index >= vec->getNumComponents()) {
		throw RuntimeError("More functions given than available components");
	}
	std::array<double, D> real_coord;
	for (int i = 0; i < vec->getNumLocalPatches(); i++) {
		auto pinfo = domain->getPatchInfoMap()[domain->patch_id_bc_map_vec[i]];
		for (Side<D> s : Side<D>::getValues()) {
			if (!pinfo.hasNbr(s)) {
				ComponentView<D - 1> ld = vec->getComponentView(component_index, pinfo.getBCLocalIndex(s));
				nested_loop<D - 1>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D - 1> &coord) {
					GetRealCoordBound<D>(pinfo, coord, s, real_coord);
					ld[coord] = func(real_coord);
				});
			}
		}
	}
}
}; // namespace DomainTools
} // namespace ThunderEgg
#endif