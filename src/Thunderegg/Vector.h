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

#ifndef THUNDEREGG_VECTOR_H
#define THUNDEREGG_VECTOR_H
#include <Thunderegg/Loops.h>
#include <Thunderegg/Side.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <numeric>
namespace Thunderegg
{
/**
 * @brief This object does any necessary cleanup when deleted
 */
class LocalDataManager
{
	public:
	virtual ~LocalDataManager(){};
};

/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template <size_t D> class LocalData
{
	private:
	double *                          data;
	std::array<int, D>                strides;
	std::array<int, D>                lengths;
	std::array<int, D>                start;
	std::array<int, D>                end;
	std::array<int, D>                ghost_start;
	std::array<int, D>                ghost_end;
	std::shared_ptr<LocalDataManager> ldm;
	int                               num_ghost_cells;

	LocalData<D - 1> getSliceOnSidePriv(Side<D> s, int offset) const;

	public:
	LocalData() {}
	/**
	 * @brief Construct a new LocalData object
	 *
	 * @param data pointer to the first element in the patch (non-ghost cell element)
	 * @param strides the strides in each direction
	 * @param lengths the lengths in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @param ldm the local data manager for the data
	 */
	LocalData(double *data, const std::array<int, D> &strides, const std::array<int, D> &lengths,
	          int num_ghost_cells, std::shared_ptr<LocalDataManager> ldm = nullptr)
	{
		this->data            = data;
		this->strides         = strides;
		this->lengths         = lengths;
		this->ldm             = ldm;
		this->num_ghost_cells = num_ghost_cells;
		start.fill(0);
		end = lengths;
		for (size_t i = 0; i < D; i++) {
			end[i]--;
		}
		ghost_start = start;
		ghost_end   = end;
		for (size_t i = 0; i < D; i++) {
			ghost_start[i] -= num_ghost_cells;
			ghost_end[i] += num_ghost_cells;
		}
	}
	/**
	 * @brief Get the pointer the data at the specified coordinate
	 *
	 * @param coord the coordianate
	 * @return double* the pointer
	 */
	inline double *getPtr(const std::array<int, D> &coord)
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data + idx;
	}
	/**
	 * @brief Get the pointer the data at the specified coordinate
	 *
	 * @param coord the coordianate
	 * @return double* the pointer
	 */
	inline const double *getPtr(const std::array<int, D> &coord) const
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data + idx;
	}
	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * @param coord the coordinate
	 * @return double& the element
	 */
	inline double &operator[](const std::array<int, D> &coord)
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data[idx];
	}
	/**
	 * @brief Get a reference to the element at the specified coordinate
	 *
	 * @param coord the coordinate
	 * @return double& the element
	 */
	inline const double &operator[](const std::array<int, D> &coord) const
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data[idx];
	}
	/**
	 * @brief Get a slice with dimensions D-1 on the specified side of the patch
	 *
	 * @param s the side
	 * @param offset how far from the side the slice is
	 * @return LocalData<D - 1>
	 */
	LocalData<D - 1> getSliceOnSide(Side<D> s, int offset = 0)
	{
		return getSliceOnSidePriv(s, offset);
	}
	/**
	 * @brief Get a slice with dimensions D-1 on the specified side of the patch
	 *
	 * @param s the side
	 * @param offset how far from the side the slice is
	 * @return LocalData<D - 1>
	 */
	const LocalData<D - 1> getSliceOnSide(Side<D> s, int offset = 0) const
	{
		return getSliceOnSidePriv(s, offset);
	}
	/**
	 * @brief Get a slice of ghost cells with dimensions D-1 on the specified side of the patch
	 *
	 * @param s the side
	 * @param offset which layer of ghost cells to acess
	 * @return LocalData<D - 1>
	 */
	LocalData<D - 1> getGhostSliceOnSide(Side<D> s, int offset) const
	{
		return getSliceOnSidePriv(s, -offset);
	}
	/**
	 * @brief Get the Lengths of the patch in each direction
	 */
	const std::array<int, D> &getLengths() const
	{
		return lengths;
	}
	/**
	 * @brief Get the strides of the patch in each direction
	 */
	const std::array<int, D> &getStrides() const
	{
		return strides;
	}
	/**
	 * @brief Get the coordinate of the first element
	 */
	const std::array<int, D> &getStart() const
	{
		return start;
	}
	/**
	 * @brief Get the coordinate of the last element
	 */
	const std::array<int, D> &getEnd() const
	{
		return end;
	}
	/**
	 * @brief Get the coordinate of the first ghost cell element
	 */
	const std::array<int, D> &getGhostStart() const
	{
		return ghost_start;
	}
	/**
	 * @brief Get the coordinate of the last ghost cell element
	 */
	const std::array<int, D> &getGhostEnd() const
	{
		return ghost_end;
	}
	/**
	 * @brief Get the number of ghost cells on each side of the patch
	 */
	int getNumGhostCells() const
	{
		return num_ghost_cells;
	}
	/**
	 * @brief Get the pointer to the first element
	 */
	double *getPtr() const
	{
		return data;
	}
};
template <size_t D>
inline LocalData<D - 1> LocalData<D>::getSliceOnSidePriv(Side<D> s, int offset) const
{
	size_t                 axis = s.toInt() / 2;
	std::array<int, D - 1> new_strides;
	for (size_t i = 0; i < axis; i++) {
		new_strides[i] = strides[i];
	}
	for (size_t i = axis; i < D - 1; i++) {
		new_strides[i] = strides[i + 1];
	}
	std::array<int, D - 1> new_lengths;
	for (size_t i = 0; i < axis; i++) {
		new_lengths[i] = lengths[i];
	}
	for (size_t i = axis; i < D - 1; i++) {
		new_lengths[i] = lengths[i + 1];
	}
	if (s.isLowerOnAxis()) {
		double *new_data = data + offset * strides[axis];
		return LocalData<D - 1>(new_data, new_strides, new_lengths, 0, ldm);
	} else {
		double *new_data = data + (lengths[axis] - 1 - offset) * strides[axis];
		return LocalData<D - 1>(new_data, new_strides, new_lengths, 0, ldm);
	}
}

/**
 * @brief Vector class for use in thunderegg
 *
 * @tparam D the number of cartesian dimensions
 */
template <size_t D> class Vector
{
	protected:
	/**
	 * @brief The number of local patches in the vector
	 */
	int num_local_patches;
	/**
	 * @brief the mpi comm
	 */
	MPI_Comm comm = MPI_COMM_WORLD;

	public:
	/**
	 * @brief Destroy the Vector object
	 */
	virtual ~Vector(){};
	/**
	 * @brief Get the number of local patches
	 */
	int getNumLocalPatches()
	{
		return num_local_patches;
	}
	/**
	 * @brief Get the LocalData object for the specified path
	 *
	 * @param patch_local_index the local index of the patch
	 * @return LocalData<D> the LocalData object
	 */
	virtual LocalData<D> getLocalData(int patch_local_index) = 0;
	/**
	 * @brief Get the LocalData object for the specified path
	 *
	 * @param patch_local_index the local index of the patch
	 * @return LocalData<D> the LocalData object
	 */
	virtual const LocalData<D> getLocalData(int patch_local_index) const = 0;
	virtual void               setNumGhostPatches(int num_ghost_patches) = 0;

	/**
	 * @brief set all value in the vector
	 *
	 * @param alpha the value ot be set
	 */
	virtual void set(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] = alpha; });
		}
	}
	/**
	 * @brief scale all elements in the vector
	 *
	 * @param alpha the value to scale by
	 */
	virtual void scale(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] *= alpha; });
		}
	}
	/**
	 * @brief shift all the values in the vector
	 *
	 * @param delta the value to shift by
	 */
	virtual void shift(double delta)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] += delta; });
		}
	}
	/**
	 * @brief copy the values of the other vector
	 *
	 * @param b the other vector
	 */
	virtual void copy(std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] = ld_b[coord]; });
		}
	}
	/**
	 * @brief add the other vector to this vector
	 *
	 * @param b the other vector
	 */
	virtual void add(std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] += ld_b[coord]; });
		}
	}
	/**
	 * @brief `this = this + alpha * b`
	 */
	virtual void addScaled(double alpha, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] += ld_b[coord] * alpha; });
		}
	}
	/**
	 * @brief `this = this + alpha * a + beta * b`
	 */
	virtual void addScaled(double alpha, std::shared_ptr<const Vector<D>> a, double beta,
	                       std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_a = a->getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) {
				ld[coord] += ld_a[coord] * alpha + ld_b[coord] * beta;
			});
		}
	}
	/**
	 * @brief `this = alpha * this + b`
	 */
	virtual void scaleThenAdd(double alpha, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) {
				ld[coord] = alpha * ld[coord] + ld_b[coord];
			});
		}
	}
	/**
	 * @brief `this = alpha * this + beta * b`
	 */
	virtual void scaleThenAddScaled(double alpha, double beta, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) {
				ld[coord] = alpha * ld[coord] + beta * ld_b[coord];
			});
		}
	}
	/**
	 * @brief `this = alpha * this + beta * b + gamma * c`
	 */
	virtual void scaleThenAddScaled(double alpha, double beta, std::shared_ptr<const Vector<D>> b,
	                                double gamma, std::shared_ptr<const Vector<D>> c)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			const LocalData<D> ld_c = c->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) {
				ld[coord] = alpha * ld[coord] + beta * ld_b[coord] + gamma * ld_c[coord];
			});
		}
	}
	/**
	 * @brief get the l2norm
	 */
	virtual double twoNorm() const
	{
		double sum = 0;
		for (int i = 0; i < num_local_patches; i++) {
			const LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { sum += ld[coord] * ld[coord]; });
		}
		double global_sum;
		MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
		return sqrt(global_sum);
	}
	/**
	 * @brief get the infnorm
	 */
	virtual double infNorm() const
	{
		double max = 0;
		for (int i = 0; i < num_local_patches; i++) {
			const LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { max = fmax(fabs(ld[coord]), max); });
		}
		double global_max;
		MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
		return global_max;
	}
	/**
	 * @brief get the infnorm
	 */
	virtual double dot(std::shared_ptr<const Vector<D>> b) const
	{
		double retval = 0;
		for (int i = 0; i < num_local_patches; i++) {
			const LocalData<D> ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { retval += ld[coord] * ld_b[coord]; });
		}
		double global_retval;
		MPI_Allreduce(&retval, &global_retval, 1, MPI_DOUBLE, MPI_SUM, comm);
		return global_retval;
	}
};
template <size_t D> class VectorGenerator
{
	public:
	virtual std::shared_ptr<Vector<D>> getNewVector() = 0;
};
extern template class LocalData<1>;
extern template class LocalData<2>;
extern template class LocalData<3>;
extern template class Vector<1>;
extern template class Vector<2>;
extern template class Vector<3>;
} // namespace Thunderegg
#endif
