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

#include <ThunderEgg/PatchInfo.h>

namespace ThunderEgg {

template<int D>
requires is_supported_dimension<D>
PatchInfo<D>::PatchInfo()
{
  starts.fill(0);
  ns.fill(1);
  spacings.fill(1);
  child_ids.fill(-1);
  child_ranks.fill(-1);
}

template<int D>
requires is_supported_dimension<D>
PatchInfo<D>::PatchInfo(const PatchInfo<D>& other_pinfo)
  : id(other_pinfo.id)
  , local_index(other_pinfo.local_index)
  , global_index(other_pinfo.global_index)
  , refine_level(other_pinfo.refine_level)
  , parent_id(other_pinfo.parent_id)
  , parent_rank(other_pinfo.parent_rank)
  , child_ids(other_pinfo.child_ids)
  , child_ranks(other_pinfo.child_ranks)
  , num_ghost_cells(other_pinfo.num_ghost_cells)
  , rank(other_pinfo.rank)
  , orth_on_parent(other_pinfo.orth_on_parent)
  , ns(other_pinfo.ns)
  , starts(other_pinfo.starts)
  , spacings(other_pinfo.spacings)
{
  for (size_t i = 0; i < other_pinfo.nbr_infos.size(); i++) {
    if (other_pinfo.nbr_infos[i] != nullptr) {
      nbr_infos[i] = other_pinfo.nbr_infos[i]->clone();
    }
  }
}

template<int D>
requires is_supported_dimension<D> PatchInfo<D>
&PatchInfo<D>::operator=(const PatchInfo<D>& other_pinfo)
{
  id = other_pinfo.id;
  local_index = other_pinfo.local_index;
  global_index = other_pinfo.global_index;
  refine_level = other_pinfo.refine_level;
  parent_id = other_pinfo.parent_id;
  parent_rank = other_pinfo.parent_rank;
  child_ids = other_pinfo.child_ids;
  child_ranks = other_pinfo.child_ranks;
  num_ghost_cells = other_pinfo.num_ghost_cells;
  rank = other_pinfo.rank;
  orth_on_parent = other_pinfo.orth_on_parent;
  ns = other_pinfo.ns;
  starts = other_pinfo.starts;
  spacings = other_pinfo.spacings;

  for (size_t i = 0; i < other_pinfo.nbr_infos.size(); i++) {
    if (other_pinfo.nbr_infos[i] != nullptr) {
      nbr_infos[i] = other_pinfo.nbr_infos[i]->clone();
    } else {
      nbr_infos[i] = nullptr;
    }
  }
  return *this;
}

template<int D>
requires is_supported_dimension<D>
bool
operator<(const PatchInfo<D>& l, const PatchInfo<D>& r)
{
  return l.id < r.id;
}

template<int D>
requires is_supported_dimension<D>
template<int M>
void
PatchInfo<D>::setNbrInfo(Face<D, M> f, std::nullptr_t)
{
  nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()] = nullptr;
}

template<int D>
requires is_supported_dimension<D>
template<int M>
void
PatchInfo<D>::setNbrInfo(Face<D, M> f, NbrInfo<M>* nbr_info)
{
  nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()].reset(nbr_info);
}

template class PatchInfo<2>;

template void
PatchInfo<2>::setNbrInfo(Face<2, 1> f, std::nullptr_t);
template void
PatchInfo<2>::setNbrInfo(Face<2, 0> f, std::nullptr_t);

template void
PatchInfo<2>::setNbrInfo(Face<2, 1> f, NbrInfo<1>* nbr_info);

template class PatchInfo<3>;

template void
PatchInfo<3>::setNbrInfo(Face<3, 2> f, std::nullptr_t);
template void
PatchInfo<3>::setNbrInfo(Face<3, 1> f, std::nullptr_t);
template void
PatchInfo<3>::setNbrInfo(Face<3, 0> f, std::nullptr_t);

}
