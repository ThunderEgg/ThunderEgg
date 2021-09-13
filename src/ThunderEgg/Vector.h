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
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/Face.h>
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/PatchView.h>
#include <cmath>
#include <mpi.h>
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
	 * @brief the communicator
	 */
	Communicator comm;

	/**
	 * @brief the pointers to patch starts
	 */
	std::vector<double *> patch_starts;

	/**
	 * @brief the strides for each axis of the patch
	 */
	std::array<int, D + 1> strides;
	/**
	 * @brief the number of non-ghost cells in each direction of the patch
	 */
	std::array<int, D + 1> lengths;

	/**
	 * @brief The number of ghost cells
	 */
	int num_ghost_cells = 0;

	/**
	 * @brief allocated data, empty of data is not managed
	 */
	std::vector<double> data;

	/**
	 * @brief The number of local cells in the vector
	 * This exclude ghost cells
	 */
	int num_local_cells = 0;

	/**
	 * @brief determine the strides from the lengths
	 */
	void determineStrides()
	{
		int curr_stride = 1;
		for (int i = 0; i < D; i++) {
			strides[i] = curr_stride;
			curr_stride *= lengths[i] + 2 * num_ghost_cells;
		}
		strides[D] = curr_stride;
	}

	/**
	 * @brief allocate the data vector and set patch_starts
	 *
	 * @param num_local_patches number of local patches
	 */
	void allocateData(int num_local_patches)
	{
		int patch_stride = strides[D] * lengths[D];
		data.resize(patch_stride * num_local_patches);
		patch_starts.resize(num_local_patches);
		for (int i = 0; i < num_local_patches; i++) {
			patch_starts[i] = data.data() + i * patch_stride;
		}
	}

	public:
	/**
	 * @brief Construct a new Vector object of size 0
	 */
	Vector()
	{
		strides.fill(0);
		lengths.fill(0);
	}
	/**
	 * @brief Construct a new Vector object with managed memory
	 *
	 * @param comm the MPI comm that is being used
	 * @param num_components the number of components for each patch
	 * @param num_local_patches the number of local patches in this vector
	 * @param num_local_cells the number of local (non-ghost) cells in this vector
	 */
	Vector(Communicator comm, const std::array<int, D> &ns, int num_components, int num_local_patches, int num_ghost_cells)
	: comm(comm),
	  num_ghost_cells(num_ghost_cells)
	{
		for (int i = 0; i < D; i++) {
			lengths[i] = ns[i];
		}
		lengths[D]      = num_components;
		int size        = 1;
		num_local_cells = 1;
		for (int i = 0; i < D; i++) {
			strides[i] = size;
			size *= lengths[i] + 2 * num_ghost_cells;
			num_local_cells *= lengths[i];
		}
		strides[D] = size;
		size *= lengths[D];
		int patch_stride = size;
		size *= num_local_patches;
		num_local_cells *= num_local_patches;
		data.resize(size);
		patch_starts.resize(num_local_patches);
		for (int i = 0; i < num_local_patches; i++) {
			patch_starts[i] = data.data() + i * patch_stride;
		}
	}
	/**
	 * @brief Construct a new Vector object for a given domain
	 *
	 * @param domain the domain
	 * @param num_components  the number of components for each patch
	 */
	Vector(const Domain<D> &domain, int num_components)
	: comm(domain.getCommunicator()),
	  num_ghost_cells(domain.getNumGhostCells()),
	  num_local_cells(domain.getNumLocalCells())
	{
		const std::array<int, D> &ns                = domain.getNs();
		int                       num_local_patches = domain.getNumLocalPatches();
		for (int i = 0; i < D; i++) {
			lengths[i] = ns[i];
		}
		lengths[D]      = num_components;
		int size        = 1;
		num_local_cells = 1;
		for (int i = 0; i < D; i++) {
			strides[i] = size;
			size *= lengths[i] + 2 * num_ghost_cells;
			num_local_cells *= lengths[i];
		}
		strides[D] = size;
		size *= lengths[D];
		int patch_stride = size;
		size *= num_local_patches;
		num_local_cells *= num_local_patches;
		data.resize(size);
		patch_starts.resize(num_local_patches);
		for (int i = 0; i < num_local_patches; i++) {
			patch_starts[i] = data.data() + i * patch_stride;
		}
	}
	/**
	 * @brief Construct a new Vector object with unmanaged memory
	 *
	 * @param comm  the communicator
	 * @param patch_starts pointers to the starts of each patch
	 * @param strides the strides
	 * @param lengths the lengths
	 * @param num_ghost_cells  the number of ghost cells
	 */
	Vector(Communicator                  comm,
	       const std::vector<double *> & patch_starts,
	       const std::array<int, D + 1> &strides,
	       const std::array<int, D + 1> &lengths,
	       int                           num_ghost_cells)
	: comm(comm),
	  patch_starts(patch_starts),
	  strides(strides),
	  lengths(lengths),
	  num_ghost_cells(num_ghost_cells)
	{
		num_local_cells = 1;
		for (int i = 0; i < D; i++) {
			num_local_cells *= lengths[i];
		}
		num_local_cells *= patch_starts.size();
	}
	/**
	 * @brief Copy constructor
	 *
	 * will copy all values
	 *
	 * @param other the vector to copy
	 */
	Vector(const Vector<D> &other)
	: comm(other.comm),
	  lengths(other.lengths),
	  num_ghost_cells(other.num_ghost_cells),
	  num_local_cells(other.num_local_cells)
	{
		if (other.data.empty()) {
			determineStrides();
			allocateData(other.getNumLocalPatches());
			copyWithGhost(other);
		} else {
			strides          = other.strides;
			data             = other.data;
			int patch_stride = data.size() / other.getNumLocalPatches();
			patch_starts.resize(other.patch_starts.size());
			for (int i = 0; i < patch_starts.size(); i++) {
				patch_starts[i] = data.data() + i * patch_stride;
			}
		}
	}
	/**
	 * @brief Copy assignment
	 *
	 * will copy all values
	 *
	 * @param other the vector to copy
	 * @return Vector<D>& this
	 */
	Vector<D> &operator=(const Vector<D> &other)
	{
		comm            = other.comm;
		lengths         = other.lengths;
		num_ghost_cells = other.num_ghost_cells;
		num_local_cells = other.num_local_cells;
		if (other.data.empty()) {
			determineStrides();
			allocateData(other.getNumLocalPatches());
			copyWithGhost(other);
		} else {
			strides          = other.strides;
			data             = other.data;
			int patch_stride = data.size() / other.getNumLocalPatches();
			patch_starts.resize(other.patch_starts.size());
			for (int i = 0; i < patch_starts.size(); i++) {
				patch_starts[i] = data.data() + i * patch_stride;
			}
		}
		return *this;
	}
	/**
	 * @brief Move constructor
	 *
	 * @param other the vector to move
	 */
	Vector(Vector<D> &&other)
	: comm(std::exchange(other.comm, Communicator())),
	  patch_starts(std::exchange(other.patch_starts, std::vector<double *>())),
	  num_ghost_cells(std::exchange(other.num_ghost_cells, 0)),
	  data(std::exchange(other.data, std::vector<double>())),
	  num_local_cells(std::exchange(other.num_local_cells, 0))
	{
		lengths.fill(0);
		std::swap(lengths, other.lengths);
		strides.fill(0);
		std::swap(strides, other.strides);
	}
	/**
	 * @brief Move assignment
	 *
	 * @param other the vector to move
	 * @return Vector<D>&& this
	 */
	Vector<D> &operator=(Vector<D> &&other)
	{
		std::swap(comm, other.comm);
		std::swap(lengths, other.lengths);
		std::swap(num_ghost_cells, other.num_ghost_cells);
		std::swap(num_local_cells, other.num_local_cells);
		std::swap(strides, other.strides);
		std::swap(data, other.data);
		std::swap(patch_starts, other.patch_starts);
		return *this;
	}
	/**
	 * @brief get the MPI Comm that this vector uses
	 *
	 * @return MPI_Comm the comm
	 */
	const Communicator &getCommunicator() const
	{
		return comm;
	}
	int getNumComponents() const
	{
		return lengths[D];
	}
	/**
	 * @brief Get the number of local patches
	 */
	int getNumLocalPatches() const
	{
		return patch_starts.size();
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
	 * @brief Get the number of ghost cells
	 *
	 * @return int the number of ghost cells
	 */
	int getNumGhostCells() const
	{
		return num_ghost_cells;
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
	PatchView<double, D> getPatchView(int patch_local_index)
	{
		if constexpr (ENABLE_DEBUG) {
			if (patch_local_index < 0 || patch_local_index >= getNumLocalPatches()) {
				throw RuntimeError("invalid patch index");
			}
		}
		return PatchView<double, D>(patch_starts[patch_local_index], strides, lengths, num_ghost_cells);
	}

	/**
	 * @brief Get the View objects for the specified patch
	 * index of View object will correspond to component index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return View<D> the View object
	 */
	PatchView<const double, D> getPatchView(int patch_local_index) const
	{
		if constexpr (ENABLE_DEBUG) {
			if (patch_local_index < 0 || patch_local_index >= getNumLocalPatches()) {
				throw RuntimeError("invalid patch index");
			}
		}
		return PatchView<const double, D>(patch_starts[patch_local_index], strides, lengths, num_ghost_cells);
	}

	/**
	 * @brief set all value in the vector
	 *
	 * @param alpha the value ot be set
	 */
	void set(double alpha)
	{
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] = b_view[coord]; });
		}
	}
	/**
	 * @brief copy the values of the other vector include ghost cell values
	 *
	 * @param b the other vector
	 */
	void copyWithGhost(const Vector<D> &b)
	{
		for (int i = 0; i < getNumLocalPatches(); i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_all_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] = b_view[coord]; });
		}
	}
	/**
	 * @brief add the other vector to this vector
	 *
	 * @param b the other vector
	 */
	void add(const Vector<D> &b)
	{
		for (int i = 0; i < getNumLocalPatches(); i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] += b_view[coord]; });
		}
	}
	/**
	 * @brief `this = this + alpha * b`
	 */
	void addScaled(double alpha, const Vector<D> &b)
	{
		for (int i = 0; i < getNumLocalPatches(); i++) {
			PatchView<double, D>       view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { view[coord] += b_view[coord] * alpha; });
		}
	}
	/**
	 * @brief `this = this + alpha * a + beta * b`
	 */
	void addScaled(double alpha, const Vector<D> &a, double beta, const Vector<D> &b)
	{
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
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
		for (int i = 0; i < getNumLocalPatches(); i++) {
			PatchView<const double, D> view = getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { sum += view[coord] * view[coord]; });
		}
		double global_sum;
		MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm.getMPIComm());
		return sqrt(global_sum);
	}
	/**
	 * @brief get the infnorm
	 */
	double infNorm() const
	{
		double max = 0;
		for (int i = 0; i < getNumLocalPatches(); i++) {
			PatchView<const double, D> view = getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { max = fmax(view[coord], max); });
		}
		double global_max;
		MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm.getMPIComm());
		return global_max;
	}
	/**
	 * @brief get the dot product
	 */
	double dot(const Vector<D> &b) const
	{
		double retval = 0;
		for (int i = 0; i < getNumLocalPatches(); i++) {
			PatchView<const double, D> view   = getPatchView(i);
			PatchView<const double, D> b_view = b.getPatchView(i);
			loop_over_interior_indexes<D + 1>(view, [&](const std::array<int, D + 1> &coord) { retval += view[coord] * b_view[coord]; });
		}
		double global_retval;
		MPI_Allreduce(&retval, &global_retval, 1, MPI_DOUBLE, MPI_SUM, comm.getMPIComm());
		return global_retval;
	}
	/**
	 * @brief Get a vector of the same length initialized to zero
	 *
	 * @return Vector<D> the vector of the same length initialize to zero
	 */
	Vector<D> getZeroClone() const
	{
		Vector<D> clone;
		clone.comm            = comm;
		clone.lengths         = lengths;
		clone.num_ghost_cells = num_ghost_cells;
		clone.num_local_cells = num_local_cells;
		clone.determineStrides();
		clone.allocateData(getNumLocalPatches());
		return clone;
	}
};
extern template class Vector<1>;
extern template class Vector<2>;
extern template class Vector<3>;
} // namespace ThunderEgg
#endif
