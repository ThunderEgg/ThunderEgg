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

#ifndef THUNDEREGG_VECTOR_H
#define THUNDEREGG_VECTOR_H
#include <ThunderEgg/Face.h>
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/PatchView.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <numeric>
namespace ThunderEgg
{
/**
 * @brief Vector class for use in thunderegg
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D> class Vector
{
	private:
	/**
	 * @brief the pointers to patch starts
	 */
	std::vector<double *> patch_starts;

	/**
	 * @brief the number of non-ghost cells in each direction of the patch
	 */
	std::array<int, D + 1> lengths;
	/**
	 * @brief the strides for each axis of the patch
	 */
	std::array<int, D + 1> strides;

	double *data = nullptr;

	/**
	 * @brief The number of ghost cells
	 */
	int num_ghost_cells;

	/**
	 * @brief The number of components for a patch
	 */
	int num_components;
	/**
	 * @brief The number of local patches in the vector
	 */
	int num_local_patches;
	/**
	 * @brief The number of local cells in the vector
	 * This exclude ghost cells
	 */
	int num_local_cells;
	/**
	 * @brief the mpi comm
	 */
	MPI_Comm comm;

	public:
	/**
	 * @brief Construct a new Vector object
	 *
	 * @param comm the MPI comm that is being used
	 * @param num_components the number of components for each patch
	 * @param num_local_patches the number of local patches in this vector
	 * @param num_local_cells the number of local (non-ghost) cells in this vector
	 */
	explicit Vector(MPI_Comm comm, int num_components, int num_local_patches, int num_local_cells)
	: num_components(num_components),
	  num_local_patches(num_local_patches),
	  num_local_cells(num_local_cells),
	  comm(comm)
	{
	}
	/*
	Vector(MPI_Comm comm, const Domain<D> &domain, int num_components)
	{
	    //
	}
	Vector(MPI_Comm comm, const Domain<D> &domain, int num_components)
	{
	    //
	}
	*/
	/**
	 * @brief Destroy the Vector object
	 */
	virtual ~Vector(){};
	/**
	 * @brief get the MPI Comm that this vector uses
	 *
	 * @return MPI_Comm the comm
	 */
	MPI_Comm getMPIComm() const
	{
		return comm;
	}
	int getNumComponents() const
	{
		return num_components;
	}
	/**
	 * @brief Get the number of local patches
	 */
	int getNumLocalPatches() const
	{
		return num_local_patches;
	}
	/**
	 * @brief Get the number of local cells int he vector (excluding ghost cells)
	 *
	 * @return int the number of local cells
	 */
	int getNumLocalCells() const
	{
		return num_local_cells;
	}
	/**
	 * @brief Get the ComponentView for the specified patch and component
	 *
	 * @param component_index the index of the component access
	 * @param patch_local_index the local index of the patch
	 * @return ComponentView<D> the View object
	 */
	ComponentView<double, D> getComponentView(int component_index, int patch_local_index)
	{
		return getPatchView(patch_local_index).getComponentView(component_index);
	}
	/**
	 * @brief Get the ComponentView for the specified patch and component
	 *
	 * @param component_index the index of the component access
	 * @param patch_local_index the local index of the patch
	 * @return ComponentView<D> the View object
	 */
	ComponentView<const double, D> getComponentView(int component_index, int patch_local_index) const
	{
		return getPatchView(patch_local_index).getComponentView(component_index);
	}
	/**
	 * @brief Get the View objects for the specified patch
	 * index of View object will correspond to component index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return View<D> the View object
	 */
	virtual PatchView<double, D> getPatchView(int patch_local_index) = 0;

	/**
	 * @brief Get the View objects for the specified patch
	 * index of View object will correspond to component index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return View<D> the View object
	 */
	virtual PatchView<const double, D> getPatchView(int patch_local_index) const = 0;

	/**
	 * @brief set all value in the vector
	 *
	 * @param alpha the value ot be set
	 */
	void set(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D> view = getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] = alpha; });
		}
	}
	/**
	 * @brief set all values in the vector (including ghost cells)
	 *
	 * @param alpha the value ot be set
	 */
	void setWithGhost(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D> view = getPatchView(i);
			loop_over_all_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] = alpha; });
		}
	}
	/**
	 * @brief scale all elements in the vector
	 *
	 * @param alpha the value to scale by
	 */
	void scale(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D> view = getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] *= alpha; });
		}
	}
	/**
	 * @brief shift all the values in the vector
	 *
	 * @param delta the value to shift by
	 */
	void shift(double delta)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D> view = getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] += delta; });
		}
	}
	/**
	 * @brief copy the values of the other vector
	 *
	 * @param b the other vector
	 */
	void copy(const Vector<D> &b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] = b_view[coord]; });
		}
	}
	/**
	 * @brief add the other vector to this vector
	 *
	 * @param b the other vector
	 */
	void add(const Vector<D> &b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] += b_view[coord]; });
		}
	}
	/**
	 * @brief `this = this + alpha * b`
	 */
	virtual void addScaled(double alpha, const Vector<D> &b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] += b_view[coord] * alpha; });
		}
	}
	/**
	 * @brief `this = this + alpha * a + beta * b`
	 */
	virtual void addScaled(double alpha, const Vector<D> &a, double beta, const Vector<D> &b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> a_view = a.getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(
			view, [&](const std::array<int, D + 1> &coord) { view[coord] += a_view[coord] * alpha + b_view[coord] * beta; });
		}
	}
	/**
	 * @brief `this = alpha * this + b`
	 */
	void scaleThenAdd(double alpha, const Vector<D> &b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] = view[coord] * alpha + b_view[coord]; });
		}
	}
	/**
	 * @brief `this = alpha * this + beta * b`
	 */
	void scaleThenAddScaled(double alpha, double beta, const Vector<D> &b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view,
			                                  [&](const std::array<int, D + 1> &coord) { view[coord] = view[coord] * alpha + b_view[coord] * beta; });
		}
	}
	/**
	 * @brief `this = alpha * this + beta * b + gamma * c`
	 */
	void scaleThenAddScaled(double alpha, double beta, const Vector<D> &b, double gamma, const Vector<D> &c)
	{
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			PatchView<const double, D> c_view = c.getPatchView(i);
			loop_over_interior_indexes<D + 1>(
			view, [&](const std::array<int, D + 1> &coord) { view[coord] = view[coord] * alpha + b_view[coord] * beta + c_view[coord] * gamma; });
		}
	}
	/**
	 * @brief get the l2norm
	 */
	double twoNorm() const
	{
		double sum = 0;
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<const double, D> view = getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { sum += view[coord] * view[coord]; });
		}
		double global_sum;
		MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
		return sqrt(global_sum);
	}
	/**
	 * @brief get the infnorm
	 */
	double infNorm() const
	{
		double max = 0;
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<const double, D> view = getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { max = fmax(view[coord], max); });
		}
		double global_max;
		MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
		return global_max;
	}
	/**
	 * @brief get the dot product
	 */
	double dot(const Vector<D> &b) const
	{
		double retval = 0;
		for (int i = 0; i < num_local_patches; i++) {
			PatchView<const double, D> view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { retval += view[coord] * b_view[coord]; });
		}
		double global_retval;
		MPI_Allreduce(&retval, &global_retval, 1, MPI_DOUBLE, MPI_SUM, comm);
		return global_retval;
	}
};
extern template class Vector<1>;
extern template class Vector<2>;
extern template class Vector<3>;
} // namespace ThunderEgg
#endif
