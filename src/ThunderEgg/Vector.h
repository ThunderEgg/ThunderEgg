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
#include <ThunderEgg/ComponentView.h>
#include <ThunderEgg/Face.h>
#include <ThunderEgg/Loops.h>
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
	 * @return View<D> the View object
	 */
	virtual ComponentView<D> getComponentView(int component_index, int patch_local_index) = 0;
	/**
	 * @brief Get the ComponentView for the specified patch and component
	 *
	 * @param component_index the index of the component access
	 * @param patch_local_index the local index of the patch
	 * @return View<D> the View object
	 */
	virtual const ComponentView<D> getComponentView(int component_index, int patch_local_index) const = 0;
	/**
	 * @brief Get the View objects for the specified patch
	 * index of View object will correspond to component index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return View<D> the View object
	 */
	std::vector<ComponentView<D>> getComponentViews(int patch_local_index)
	{
		std::vector<ComponentView<D>> local_datas;
		local_datas.reserve(num_components);
		for (int c = 0; c < num_components; c++) {
			local_datas.emplace_back(std::move(getComponentView(c, patch_local_index)));
		}
		return local_datas;
	}
	/**
	 * @brief Get the View objects for the specified patch
	 * index of View object will correspond to component index
	 *
	 * @param patch_local_index the local index of the patch
	 * @return View<D> the View object
	 */
	const std::vector<ComponentView<D>> getComponentViews(int patch_local_index) const
	{
		std::vector<ComponentView<D>> local_datas;
		local_datas.reserve(num_components);
		for (int c = 0; c < num_components; c++) {
			local_datas.emplace_back(std::move(getComponentView(c, patch_local_index)));
		}
		return local_datas;
	}

	/**
	 * @brief set all value in the vector
	 *
	 * @param alpha the value ot be set
	 */
	virtual void set(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			std::vector<ComponentView<D>> lds = getComponentViews(i);
			for (auto &ld : lds) {
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) { ld[coord] = alpha; });
			}
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
			std::vector<ComponentView<D>> lds = getComponentViews(i);
			for (auto &ld : lds) {
				nested_loop<D>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, D> coord) { ld[coord] = alpha; });
			}
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
			std::vector<ComponentView<D>> lds = getComponentViews(i);
			for (auto &ld : lds) {
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) { ld[coord] *= alpha; });
			}
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
			std::vector<ComponentView<D>> lds = getComponentViews(i);
			for (auto &ld : lds) {
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) { ld[coord] += delta; });
			}
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
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			for (int c = 0; c < num_components; c++) {
				nested_loop<D>(lds[c].getStart(), lds[c].getEnd(), [&](std::array<int, D> coord) { lds[c][coord] = lds_b[c][coord]; });
			}
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
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			for (int c = 0; c < num_components; c++) {
				nested_loop<D>(lds[c].getStart(), lds[c].getEnd(), [&](std::array<int, D> coord) { lds[c][coord] += lds_b[c][coord]; });
			}
		}
	}
	/**
	 * @brief `this = this + alpha * b`
	 */
	virtual void addScaled(double alpha, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			for (int c = 0; c < num_components; c++) {
				nested_loop<D>(lds[c].getStart(), lds[c].getEnd(), [&](std::array<int, D> coord) { lds[c][coord] += lds_b[c][coord] * alpha; });
			}
		}
	}
	/**
	 * @brief `this = this + alpha * a + beta * b`
	 */
	virtual void addScaled(double alpha, std::shared_ptr<const Vector<D>> a, double beta, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_a = a->getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			for (int c = 0; c < num_components; c++) {
				nested_loop<D>(lds[c].getStart(), lds[c].getEnd(), [&](std::array<int, D> coord) {
					lds[c][coord] += lds_a[c][coord] * alpha + lds_b[c][coord] * beta;
				});
			}
		}
	}
	/**
	 * @brief `this = alpha * this + b`
	 */
	virtual void scaleThenAdd(double alpha, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			for (int c = 0; c < num_components; c++) {
				nested_loop<D>(
				lds[c].getStart(), lds[c].getEnd(), [&](std::array<int, D> coord) { lds[c][coord] = alpha * lds[c][coord] + lds_b[c][coord]; });
			}
		}
	}
	/**
	 * @brief `this = alpha * this + beta * b`
	 */
	virtual void scaleThenAddScaled(double alpha, double beta, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			for (int c = 0; c < num_components; c++) {
				nested_loop<D>(lds[c].getStart(), lds[c].getEnd(), [&](std::array<int, D> coord) {
					lds[c][coord] = alpha * lds[c][coord] + beta * lds_b[c][coord];
				});
			}
		}
	}
	/**
	 * @brief `this = alpha * this + beta * b + gamma * c`
	 */
	virtual void scaleThenAddScaled(double alpha, double beta, std::shared_ptr<const Vector<D>> b, double gamma, std::shared_ptr<const Vector<D>> c)
	{
		for (int i = 0; i < num_local_patches; i++) {
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			const std::vector<ComponentView<D>> lds_c = c->getComponentViews(i);
			for (int comp = 0; comp < num_components; comp++) {
				nested_loop<D>(lds[comp].getStart(), lds[comp].getEnd(), [&](std::array<int, D> coord) {
					lds[comp][coord] = alpha * lds[comp][coord] + beta * lds_b[comp][coord] + gamma * lds_c[comp][coord];
				});
			}
		}
	}
	/**
	 * @brief get the l2norm
	 */
	virtual double twoNorm() const
	{
		double sum = 0;
		for (int i = 0; i < num_local_patches; i++) {
			const std::vector<ComponentView<D>> lds = getComponentViews(i);
			for (const auto &ld : lds) {
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) { sum += ld[coord] * ld[coord]; });
			}
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
			std::vector<ComponentView<D>> lds = getComponentViews(i);
			for (const auto &ld : lds) {
				nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) { max = fmax(fabs(ld[coord]), max); });
			}
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
			std::vector<ComponentView<D>>       lds   = getComponentViews(i);
			const std::vector<ComponentView<D>> lds_b = b->getComponentViews(i);
			for (int c = 0; c < num_components; c++) {
				nested_loop<D>(lds[c].getStart(), lds[c].getEnd(), [&](std::array<int, D> coord) { retval += lds[c][coord] * lds_b[c][coord]; });
			}
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
