/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2020 Scott Aiton
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

#include "Vector.h"

namespace ThunderEgg {

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::determineStrides()
{
  int curr_stride = 1;
  for (int i = 0; i < D; i++) {
    strides[i] = curr_stride;
    curr_stride *= lengths[i] + 2 * num_ghost_cells;
  }
  strides[D] = curr_stride;
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::allocateData(int num_local_patches)
{
  int patch_stride = strides[D] * lengths[D];
  data.resize(patch_stride * num_local_patches);
  patch_starts.resize(num_local_patches);
  for (int i = 0; i < num_local_patches; i++) {
    patch_starts[i] = data.data() + i * patch_stride;
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>::Vector()
{
  strides.fill(0);
  lengths.fill(0);
}

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>::Vector(Communicator comm,
                  const std::array<int, D>& ns,
                  int num_components,
                  int num_local_patches,
                  int num_ghost_cells)
  : comm(comm)
  , num_ghost_cells(num_ghost_cells)
{
  for (int i = 0; i < D; i++) {
    lengths[i] = ns[i];
  }
  lengths[D] = num_components;
  int size = 1;
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

template<int D>
  requires is_supported_vector_dimension<D>
template<int DomainD>
Vector<D>::Vector(const Domain<DomainD>& domain, int num_components)
  requires is_supported_dimension<D> && (D == DomainD)
: comm(domain.getCommunicator())
, num_ghost_cells(domain.getNumGhostCells())
, num_local_cells(domain.getNumLocalCells())
{
  const std::array<int, D>& ns = domain.getNs();
  int num_local_patches = domain.getNumLocalPatches();
  for (int i = 0; i < D; i++) {
    lengths[i] = ns[i];
  }
  lengths[D] = num_components;
  int size = 1;
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

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>::Vector(Communicator comm,
                  const std::vector<double*>& patch_starts,
                  const std::array<int, D + 1>& strides,
                  const std::array<int, D + 1>& lengths,
                  int num_ghost_cells)
  : comm(comm)
  , patch_starts(patch_starts)
  , strides(strides)
  , lengths(lengths)
  , num_ghost_cells(num_ghost_cells)
{
  num_local_cells = 1;
  for (int i = 0; i < D; i++) {
    num_local_cells *= lengths[i];
  }
  num_local_cells *= patch_starts.size();
}

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>::Vector(const Vector<D>& other)
  : comm(other.comm)
  , lengths(other.lengths)
  , num_ghost_cells(other.num_ghost_cells)
  , num_local_cells(other.num_local_cells)
{
  if (other.data.empty()) {
    determineStrides();
    allocateData(other.getNumLocalPatches());
    copyWithGhost(other);
  } else {
    strides = other.strides;
    data = other.data;
    int patch_stride = data.size() / other.getNumLocalPatches();
    patch_starts.resize(other.patch_starts.size());
    for (int i = 0; i < patch_starts.size(); i++) {
      patch_starts[i] = data.data() + i * patch_stride;
    }
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>&
Vector<D>::operator=(const Vector<D>& other)
{
  comm = other.comm;
  lengths = other.lengths;
  num_ghost_cells = other.num_ghost_cells;
  num_local_cells = other.num_local_cells;
  if (other.data.empty()) {
    determineStrides();
    allocateData(other.getNumLocalPatches());
    copyWithGhost(other);
  } else {
    strides = other.strides;
    data = other.data;
    int patch_stride = data.size() / other.getNumLocalPatches();
    patch_starts.resize(other.patch_starts.size());
    for (int i = 0; i < patch_starts.size(); i++) {
      patch_starts[i] = data.data() + i * patch_stride;
    }
  }
  return *this;
}

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>::Vector(Vector<D>&& other)
  : comm(std::exchange(other.comm, Communicator()))
  , patch_starts(std::exchange(other.patch_starts, std::vector<double*>()))
  , num_ghost_cells(std::exchange(other.num_ghost_cells, 0))
  , data(std::exchange(other.data, std::vector<double>()))
  , num_local_cells(std::exchange(other.num_local_cells, 0))
{
  lengths.fill(0);
  std::swap(lengths, other.lengths);
  strides.fill(0);
  std::swap(strides, other.strides);
}

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>&
Vector<D>::operator=(Vector<D>&& other)
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

template<int D>
  requires is_supported_vector_dimension<D>
const Communicator&
Vector<D>::getCommunicator() const
{
  return comm;
}

template<int D>
  requires is_supported_vector_dimension<D>
int
Vector<D>::getNumComponents() const
{
  return lengths[D];
}

template<int D>
  requires is_supported_vector_dimension<D>
int
Vector<D>::getNumLocalPatches() const
{
  return patch_starts.size();
}

template<int D>
  requires is_supported_vector_dimension<D>
int
Vector<D>::getNumLocalCells() const
{
  return num_local_cells;
}

template<int D>
  requires is_supported_vector_dimension<D>
int
Vector<D>::getNumGhostCells() const
{
  return num_ghost_cells;
}

template<int D>
  requires is_supported_vector_dimension<D>
ComponentView<double, D>
Vector<D>::getComponentView(int component_index, int patch_local_index)
{
  return getPatchView(patch_local_index).getComponentView(component_index);
}

template<int D>
  requires is_supported_vector_dimension<D>
ComponentView<const double, D>
Vector<D>::getComponentView(int component_index, int patch_local_index) const
{
  return getPatchView(patch_local_index).getComponentView(component_index);
}

template<int D>
  requires is_supported_vector_dimension<D>
PatchView<double, D>
Vector<D>::getPatchView(int patch_local_index)
{
  if constexpr (ENABLE_DEBUG) {
    if (patch_local_index < 0 || patch_local_index >= getNumLocalPatches()) {
      throw RuntimeError("invalid patch index");
    }
  }
  return PatchView<double, D>(patch_starts[patch_local_index], strides, lengths, num_ghost_cells);
}

template<int D>
  requires is_supported_vector_dimension<D>
PatchView<const double, D>
Vector<D>::getPatchView(int patch_local_index) const
{
  if constexpr (ENABLE_DEBUG) {
    if (patch_local_index < 0 || patch_local_index >= getNumLocalPatches()) {
      throw RuntimeError("invalid patch index");
    }
  }
  return PatchView<const double, D>(
    patch_starts[patch_local_index], strides, lengths, num_ghost_cells);
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::set(double alpha)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { view[coord] = alpha; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::setWithGhost(double alpha)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    Loop::OverAllIndexes<D + 1>(view,
                                [&](const std::array<int, D + 1>& coord) { view[coord] = alpha; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::scale(double alpha)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { view[coord] *= alpha; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::shift(double delta)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { view[coord] += delta; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::copy(const Vector<D>& b)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { view[coord] = b_view[coord]; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::copyWithGhost(const Vector<D>& b)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverAllIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { view[coord] = b_view[coord]; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::add(const Vector<D>& b)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { view[coord] += b_view[coord]; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::addScaled(double alpha, const Vector<D>& b)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { view[coord] += b_view[coord] * alpha; });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::addScaled(double alpha, const Vector<D>& a, double beta, const Vector<D>& b)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> a_view = a.getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(view, [&](const std::array<int, D + 1>& coord) {
      view[coord] += a_view[coord] * alpha + b_view[coord] * beta;
    });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::scaleThenAdd(double alpha, const Vector<D>& b)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(view, [&](const std::array<int, D + 1>& coord) {
      view[coord] = view[coord] * alpha + b_view[coord];
    });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::scaleThenAddScaled(double alpha, double beta, const Vector<D>& b)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(view, [&](const std::array<int, D + 1>& coord) {
      view[coord] = view[coord] * alpha + b_view[coord] * beta;
    });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
void
Vector<D>::scaleThenAddScaled(double alpha,
                              double beta,
                              const Vector<D>& b,
                              double gamma,
                              const Vector<D>& c)
{
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    PatchView<const double, D> c_view = c.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(view, [&](const std::array<int, D + 1>& coord) {
      view[coord] = view[coord] * alpha + b_view[coord] * beta + c_view[coord] * gamma;
    });
  }
}

template<int D>
  requires is_supported_vector_dimension<D>
double
Vector<D>::twoNorm() const
{
  double sum = 0;
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<const double, D> view = getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { sum += view[coord] * view[coord]; });
  }
  double global_sum;
  MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm.getMPIComm());
  return sqrt(global_sum);
}

template<int D>
  requires is_supported_vector_dimension<D>
double
Vector<D>::infNorm() const
{
  double max = 0;
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<const double, D> view = getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { max = fmax(view[coord], max); });
  }
  double global_max;
  MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm.getMPIComm());
  return global_max;
}

template<int D>
  requires is_supported_vector_dimension<D>
double
Vector<D>::dot(const Vector<D>& b) const
{
  double retval = 0;
  for (int i = 0; i < getNumLocalPatches(); i++) {
    PatchView<const double, D> view = getPatchView(i);
    PatchView<const double, D> b_view = b.getPatchView(i);
    Loop::OverInteriorIndexes<D + 1>(
      view, [&](const std::array<int, D + 1>& coord) { retval += view[coord] * b_view[coord]; });
  }
  double global_retval;
  MPI_Allreduce(&retval, &global_retval, 1, MPI_DOUBLE, MPI_SUM, comm.getMPIComm());
  return global_retval;
}

template<int D>
  requires is_supported_vector_dimension<D>
Vector<D>
Vector<D>::getZeroClone() const
{
  Vector<D> clone;
  clone.comm = comm;
  clone.lengths = lengths;
  clone.num_ghost_cells = num_ghost_cells;
  clone.num_local_cells = num_local_cells;
  clone.determineStrides();
  clone.allocateData(getNumLocalPatches());
  return clone;
}

// EXPLICIT INSTANTIATIONS

template class Vector<1>;

template class Vector<2>;
template Vector<2>::Vector(const Domain<2>& domain, int num_components);

template class Vector<3>;
template Vector<3>::Vector(const Domain<3>& domain, int num_components);

} // namespace ThunderEgg
