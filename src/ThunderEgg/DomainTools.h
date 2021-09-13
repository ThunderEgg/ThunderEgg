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

#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Vector.h>

namespace ThunderEgg
{
/**
 * @brief Various tools for filling in values in a domain
 */
class DomainTools
{
	private:
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
	template <int D, typename T> static void _SetValues(const Domain<D> &domain, Vector<D> &vec, int component_index, T func)
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
	static void _SetValues(const Domain<D> &domain, Vector<D> &vec, int component_index, T func, Args... args)
	{
		SetValues(domain, vec, component_index, func);
		_SetValues(domain, vec, component_index + 1, args...);
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
	template <int D, typename T> static void _SetValuesWithGhost(const Domain<D> &domain, Vector<D> &vec, int component_index, T func)
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
	static void _SetValuesWithGhost(const Domain<D> &domain, Vector<D> &vec, int component_index, T func, Args... args)
	{
		SetValuesWithGhost(domain, vec, component_index, func);
		_SetValuesWithGhost(domain, vec, component_index + 1, args...);
	}

	public:
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
	template <int D> static void GetRealCoord(const PatchInfo<D> &pinfo, const std::array<int, D> &coord, std::array<double, D> &real_coord)
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
	template <int D> static void GetRealCoordGhost(const PatchInfo<D> &pinfo, const std::array<int, D> &coord, std::array<double, D> &real_coord)
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
	template <int D>
	static void GetRealCoordBound(const PatchInfo<D> &pinfo, const std::array<int, D - 1> &coord, Side<D> s, std::array<double, D> &real_coord)
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
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	template <int D>
	static void
	SetValues(const Domain<D> &domain, Vector<D> &vec, int component_index, std::function<double(const std::array<double, (int) D> &)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		std::array<double, D> real_coord;
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, D> ld    = vec.getComponentView(component_index, i);
			auto                     pinfo = domain.getPatchInfoVector()[i];
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D> &coord) {
				GetRealCoord<D>(pinfo, coord, real_coord);
				ld[coord] = func(real_coord);
			});
		}
	}
	/**
	 * @brief Set the values for a vector with the given function.
	 *
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	static void SetValues(const Domain<3> &domain, Vector<3> &vec, int component_index, std::function<double(double, double, double)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, 3> ld    = vec.getComponentView(component_index, i);
			const PatchInfo<3> &     pinfo = domain.getPatchInfoVector()[i];
			double                   dx    = pinfo.spacings[0];
			double                   dy    = pinfo.spacings[1];
			double                   dz    = pinfo.spacings[2];
			for (int zi = 0; zi < pinfo.ns[2]; zi++) {
				double z = pinfo.starts[2] + 0.5 * dz + zi * dz;
				for (int yi = 0; yi < pinfo.ns[1]; yi++) {
					double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
					for (int xi = 0; xi < pinfo.ns[0]; xi++) {
						double x       = pinfo.starts[0] + 0.5 * dx + xi * dx;
						ld(xi, yi, zi) = func(x, y, z);
					}
				}
			}
		}
	}
	/**
	 * @brief Set the values for a vector with the given function.
	 *
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	static void SetValues(const Domain<2> &domain, Vector<2> &vec, int component_index, std::function<double(double, double)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, 2> ld    = vec.getComponentView(component_index, i);
			const PatchInfo<2> &     pinfo = domain.getPatchInfoVector()[i];
			double                   dx    = pinfo.spacings[0];
			double                   dy    = pinfo.spacings[1];
			for (int yi = 0; yi < pinfo.ns[1]; yi++) {
				double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
				for (int xi = 0; xi < pinfo.ns[0]; xi++) {
					double x   = pinfo.starts[0] + 0.5 * dx + xi * dx;
					ld(xi, yi) = func(x, y);
				}
			}
		}
	}
	/**
	 * @brief Set the values for a vector with the given function.
	 *
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	static void SetValues(const Domain<1> &domain, Vector<1> &vec, int component_index, std::function<double(double)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, 1> ld    = vec.getComponentView(component_index, i);
			const PatchInfo<1> &     pinfo = domain.getPatchInfoVector()[i];
			double                   dx    = pinfo.spacings[0];
			for (int xi = 0; xi < pinfo.ns[0]; xi++) {
				double x = pinfo.starts[0] + 0.5 * dx + xi * dx;
				ld(xi)   = func(x);
			}
		}
	}
	/**
	 * @brief Set the values for a vector with the given functions
	 *
	 * @tparam D the number of cartesian dimensions
	 * @tparam Args additional functions
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 * @param args additional functions for additional components
	 */
	template <int D, typename... Args>
	static void SetValues(const Domain<D> &domain, Vector<D> &vec, std::function<double(const std::array<double, D> &)> func, Args... args)
	{
		_SetValues(domain, vec, 0, func, args...);
	}
	/**
	 * @brief Set the values (including ghost values) for a vector with the given function.
	 *
	 * @tparam D the number of cartesian dimensions
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	template <int D>
	static void
	SetValuesWithGhost(const Domain<D> &domain, Vector<D> &vec, int component_index, std::function<double(const std::array<double, (int) D> &)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		std::array<double, D> real_coord;
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, D> ld    = vec.getComponentView(component_index, i);
			auto                     pinfo = domain.getPatchInfoVector()[i];
			nested_loop<D>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, D> &coord) {
				GetRealCoordGhost<D>(pinfo, coord, real_coord);
				ld[coord] = func(real_coord);
			});
		}
	}
	/**
	 * @brief Set the values (including ghost values) for a vector with the given function.
	 *
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	static void SetValuesWithGhost(const Domain<3> &domain, Vector<3> &vec, int component_index, std::function<double(double, double, double)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, 3> ld        = vec.getComponentView(component_index, i);
			const PatchInfo<3> &     pinfo     = domain.getPatchInfoVector()[i];
			int                      num_ghost = pinfo.num_ghost_cells;
			double                   dx        = pinfo.spacings[0];
			double                   dy        = pinfo.spacings[1];
			double                   dz        = pinfo.spacings[2];
			for (int zi = -num_ghost; zi < pinfo.ns[2] + num_ghost; zi++) {
				double z = pinfo.starts[2] + 0.5 * dz + zi * dz;
				for (int yi = -num_ghost; yi < pinfo.ns[1] + num_ghost; yi++) {
					double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
					for (int xi = -num_ghost; xi < pinfo.ns[0] + num_ghost; xi++) {
						double x       = pinfo.starts[0] + 0.5 * dx + xi * dx;
						ld(xi, yi, zi) = func(x, y, z);
					}
				}
			}
		}
	}
	/**
	 * @brief Set the values (including ghost values) for a vector with the given function.
	 *
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	static void SetValuesWithGhost(const Domain<2> &domain, Vector<2> &vec, int component_index, std::function<double(double, double)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, 2> ld        = vec.getComponentView(component_index, i);
			const PatchInfo<2> &     pinfo     = domain.getPatchInfoVector()[i];
			int                      num_ghost = pinfo.num_ghost_cells;
			double                   dx        = pinfo.spacings[0];
			double                   dy        = pinfo.spacings[1];
			for (int yi = -num_ghost; yi < pinfo.ns[1] + num_ghost; yi++) {
				double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
				for (int xi = -num_ghost; xi < pinfo.ns[0] + num_ghost; xi++) {
					double x   = pinfo.starts[0] + 0.5 * dx + xi * dx;
					ld(xi, yi) = func(x, y);
				}
			}
		}
	}
	/**
	 * @brief Set the values (including ghost values) for a vector with the given function.
	 *
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 */
	static void SetValuesWithGhost(const Domain<1> &domain, Vector<1> &vec, int component_index, std::function<double(double)> func)
	{
		if (component_index >= vec.getNumComponents()) {
			throw RuntimeError("Invalid component to set");
		}
		for (int i = 0; i < vec.getNumLocalPatches(); i++) {
			ComponentView<double, 1> ld        = vec.getComponentView(component_index, i);
			const PatchInfo<1> &     pinfo     = domain.getPatchInfoVector()[i];
			int                      num_ghost = pinfo.num_ghost_cells;
			double                   dx        = pinfo.spacings[0];
			for (int xi = -num_ghost; xi < pinfo.ns[0] + num_ghost; xi++) {
				double x = pinfo.starts[0] + 0.5 * dx + xi * dx;
				ld(xi)   = func(x);
			}
		}
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
	template <int D, typename... Args>
	static void SetValuesWithGhost(const Domain<D> &domain, Vector<D> &vec, std::function<double(const std::array<double, D> &)> func, Args... args)
	{
		_SetValuesWithGhost(domain, vec, 0, func, args...);
	}
	/**
	 * @brief Set the values (including ghost values) for a vector with the given functions
	 *
	 * @tparam Args additional functions
	 * @param domain the domain that we are setting values for
	 * @param vec the vector to set values in
	 * @param component_index the component to set
	 * @param func the function
	 * @param args additional functions for additional components
	 */
	template <typename... Args>
	static void SetValuesWithGhost(const Domain<3> &domain, Vector<3> &vec, std::function<double(double, double, double)> func, Args... args)
	{
		_SetValuesWithGhost(domain, vec, 0, func, args...);
	}
	/**
	 * @brief Integrate a vector over the domain.
	 *
	 * @param u the vector
	 * @return double the result of the integral
	 */
	template <int D> static double Integrate(const Domain<D> &domain, const Vector<D> &u)
	{
		double sum = 0;

		for (const auto &pinfo : domain.getPatchInfoVector()) {
			for (int c = 0; c < u.getNumComponents(); c++) {
				ComponentView<const double, D> u_data = u.getComponentView(c, pinfo.local_index);

				double patch_sum = 0;
				nested_loop<D>(u_data.getStart(), u_data.getEnd(), [&](std::array<int, D> coord) { patch_sum += u_data[coord]; });

				for (size_t i = 0; i < D; i++) {
					patch_sum *= pinfo.spacings[i];
				}
				sum += patch_sum;
			}
		}
		double retval;
		MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return retval;
	}
};
} // namespace ThunderEgg
#endif