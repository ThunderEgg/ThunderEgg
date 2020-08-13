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
#include <ThunderEgg/LocalData.h>
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/Side.h>
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
template <size_t D> class Vector
{
	private:
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
	 * @param comm_in the MPI comm that is being used
	 * @param num_local_patches_in the number of local patches in this vector
	 * @param num_local_cells_in the number of local (non-ghost) cells in this vector
	 */
	explicit Vector(MPI_Comm comm_in, int num_local_patches_in, int num_local_cells_in)
	: num_local_patches(num_local_patches_in), num_local_cells(num_local_cells_in), comm(comm_in){};
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
	 * @brief set all values in the vector (including ghost cells)
	 *
	 * @param alpha the value ot be set
	 */
	virtual void setWithGhost(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getGhostStart(), ld.getGhostEnd(),
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
	 * @brief get the dot product
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
extern template class Vector<1>;
extern template class Vector<2>;
extern template class Vector<3>;
} // namespace ThunderEgg
#endif
