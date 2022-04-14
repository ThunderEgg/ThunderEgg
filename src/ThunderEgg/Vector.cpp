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
template<int DomainD>
Vector<D>::Vector(const Domain<DomainD>& domain,
                  int num_components) requires is_supported_dimension<D> &&(D == DomainD)
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

template class Vector<1>;

template class Vector<2>;
template Vector<2>::Vector(const Domain<2>& domain, int num_components);

template class Vector<3>;
template Vector<3>::Vector(const Domain<3>& domain, int num_components);

} // namespace ThunderEgg
