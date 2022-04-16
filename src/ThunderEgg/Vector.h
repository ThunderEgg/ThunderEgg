/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2021 Scott Aiton
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
/**
 * @file
 *
 * @brief Vector class
 */
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/Face.h>
#include <ThunderEgg/Loops.h>
#include <ThunderEgg/PatchView.h>
#include <cmath>
#include <mpi.h>
#include <utility>

namespace ThunderEgg {

template<int D>
constexpr bool is_supported_vector_dimension = D > 0 && D <= 3;

/**
 * @brief Vector class for use in thunderegg
 *
 * @tparam D the number of cartesian dimensions
 */
template<int D>
  requires is_supported_vector_dimension<D>
class Vector
{
private:
  /**
   * @brief the communicator
   */
  Communicator comm;

  /**
   * @brief the pointers to patch starts
   */
  std::vector<double*> patch_starts;

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
  void
  determineStrides();

  /**
   * @brief allocate the data vector and set patch_starts
   *
   * @param num_local_patches number of local patches
   */
  void
  allocateData(int num_local_patches);

public:
  /**
   * @brief Construct a new Vector object of size 0
   */
  Vector();

  /**
   * @brief Construct a new Vector object with managed memory
   *
   * @param comm the MPI comm that is being used
   * @param num_components the number of components for each patch
   * @param num_local_patches the number of local patches in this vector
   * @param num_local_cells the number of local (non-ghost) cells in this vector
   */
  Vector(Communicator comm,
         const std::array<int, D>& ns,
         int num_components,
         int num_local_patches,
         int num_ghost_cells);

  /**
   * @brief Construct a new Vector object for a given domain
   *
   * @param domain the domain
   * @param num_components  the number of components for each patch
   */
  template<int DomainD>
  Vector(const Domain<DomainD>& domain, int num_components)
    requires is_supported_dimension<D> && (D == DomainD);

  /**
   * @brief Construct a new Vector object with unmanaged memory
   *
   * @param comm  the communicator
   * @param patch_starts pointers to the starts of each patch
   * @param strides the strides
   * @param lengths the lengths
   * @param num_ghost_cells  the number of ghost cells
   */
  Vector(Communicator comm,
         const std::vector<double*>& patch_starts,
         const std::array<int, D + 1>& strides,
         const std::array<int, D + 1>& lengths,
         int num_ghost_cells);

  /**
   * @brief Copy constructor
   *
   * will copy all values
   *
   * @param other the vector to copy
   */
  Vector(const Vector<D>& other);

  /**
   * @brief Copy assignment
   *
   * will copy all values
   *
   * @param other the vector to copy
   * @return Vector<D>& this
   */
  Vector<D>&
  operator=(const Vector<D>& other);

  /**
   * @brief Move constructor
   *
   * @param other the vector to move
   */
  Vector(Vector<D>&& other);

  /**
   * @brief Move assignment
   *
   * @param other the vector to move
   * @return Vector<D>&& this
   */
  Vector<D>&
  operator=(Vector<D>&& other);

  /**
   * @brief get the MPI Comm that this vector uses
   *
   * @return MPI_Comm the comm
   */
  const Communicator&
  getCommunicator() const;

  /**
   * @brief Get the number of components
   *
   * @return int the number of components
   */
  int
  getNumComponents() const;

  /**
   * @brief Get the number of local patches
   */
  int
  getNumLocalPatches() const;

  /**
   * @brief Get the number of local cells int he vector (excluding ghost cells)
   *
   * @return int the number of local cells
   */
  int
  getNumLocalCells() const;

  /**
   * @brief Get the number of ghost cells
   *
   * @return int the number of ghost cells
   */
  int
  getNumGhostCells() const;

  /**
   * @brief Get the ComponentView for the specified patch and component
   *
   * @param component_index the index of the component access
   * @param patch_local_index the local index of the patch
   * @return ComponentView<D> the View object
   */
  ComponentView<double, D>
  getComponentView(int component_index, int patch_local_index);

  /**
   * @brief Get the ComponentView for the specified patch and component
   *
   * @param component_index the index of the component access
   * @param patch_local_index the local index of the patch
   * @return ComponentView<D> the View object
   */
  ComponentView<const double, D>
  getComponentView(int component_index, int patch_local_index) const;

  /**
   * @brief Get the View objects for the specified patch
   * index of View object will correspond to component index
   *
   * @param patch_local_index the local index of the patch
   * @return View<D> the View object
   */
  PatchView<double, D>
  getPatchView(int patch_local_index);

  /**
   * @brief Get the View objects for the specified patch
   * index of View object will correspond to component index
   *
   * @param patch_local_index the local index of the patch
   * @return View<D> the View object
   */
  PatchView<const double, D>
  getPatchView(int patch_local_index) const;

  /**
   * @brief set all value in the vector
   *
   * @param alpha the value ot be set
   */
  void
  set(double alpha);

  /**
   * @brief set all values in the vector (including ghost cells)
   *
   * @param alpha the value ot be set
   */
  void
  setWithGhost(double alpha);

  /**
   * @brief scale all elements in the vector
   *
   * @param alpha the value to scale by
   */
  void
  scale(double alpha);

  /**
   * @brief shift all the values in the vector
   *
   * @param delta the value to shift by
   */
  void
  shift(double delta);

  /**
   * @brief copy the values of the other vector
   *
   * @param b the other vector
   */
  void
  copy(const Vector<D>& b);

  /**
   * @brief copy the values of the other vector include ghost cell values
   *
   * @param b the other vector
   */
  void
  copyWithGhost(const Vector<D>& b);

  /**
   * @brief add the other vector to this vector
   *
   * @param b the other vector
   */
  void
  add(const Vector<D>& b);

  /**
   * @brief `this = this + alpha * b`
   */
  void
  addScaled(double alpha, const Vector<D>& b);

  /**
   * @brief `this = this + alpha * a + beta * b`
   */
  void
  addScaled(double alpha, const Vector<D>& a, double beta, const Vector<D>& b);

  /**
   * @brief `this = alpha * this + b`
   */
  void
  scaleThenAdd(double alpha, const Vector<D>& b);

  /**
   * @brief `this = alpha * this + beta * b`
   */
  void
  scaleThenAddScaled(double alpha, double beta, const Vector<D>& b);

  /**
   * @brief `this = alpha * this + beta * b + gamma * c`
   */
  void
  scaleThenAddScaled(double alpha,
                     double beta,
                     const Vector<D>& b,
                     double gamma,
                     const Vector<D>& c);

  /**
   * @brief get the l2norm
   */
  double
  twoNorm() const;

  /**
   * @brief get the infnorm
   */
  double
  infNorm() const;

  /**
   * @brief get the dot product
   */
  double
  dot(const Vector<D>& b) const;

  /**
   * @brief Get a vector of the same length initialized to zero
   *
   * @return Vector<D> the vector of the same length initialize to zero
   */
  Vector<D>
  getZeroClone() const;
};

// EXPLICIT INSTANTIATIONS

extern template class Vector<1>;
extern template class Vector<2>;
extern template class Vector<3>;

} // namespace ThunderEgg
#endif
