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
  requires is_supported_dimension<D>
PatchInfo<D>&
PatchInfo<D>::operator=(const PatchInfo<D>& other_pinfo)
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
  requires is_supported_dimension<D> bool
operator<(const PatchInfo<D>& l, const PatchInfo<D>& r)
{
  return l.id < r.id;
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
  requires is_valid_face<D, M>
void
PatchInfo<D>::setNbrInfo(Face<D, M> f, std::nullptr_t)
{
  nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()] = nullptr;
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
  requires is_valid_face<D, M>
void
PatchInfo<D>::setNbrInfo(Face<D, M> f, NbrInfo<M>* nbr_info)
{
  nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()].reset(nbr_info);
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
  requires is_valid_face<D, M>
NbrType
PatchInfo<D>::getNbrType(Face<D, M> s) const
{
  const NbrInfoBase* info = nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get();
  if (info == nullptr) {
    throw RuntimeError("PatchInfo::getNbrType: nbr_info is nullptr");
  }
  return info->getNbrType();
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
  requires is_valid_face<D, M>
const NormalNbrInfo<M>&
PatchInfo<D>::getNormalNbrInfo(Face<D, M> s) const
{
  const NbrInfoBase* info = nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get();
  if (info == nullptr) {
    throw RuntimeError("PatchInfo::getNormalNbrInfo: nbr_info is nullptr");
  }
  return dynamic_cast<const NormalNbrInfo<M>&>(*info);
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
  requires is_valid_face<D, M>
const CoarseNbrInfo<M>&
PatchInfo<D>::getCoarseNbrInfo(Face<D, M> s) const
{
  const NbrInfoBase* info = nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get();
  if (info == nullptr) {
    throw RuntimeError("PatchInfo::getCoarseNbrInfo: nbr_info is nullptr");
  }
  return dynamic_cast<const CoarseNbrInfo<M>&>(*info);
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
  requires is_valid_face<D, M>
const FineNbrInfo<M>&
PatchInfo<D>::getFineNbrInfo(Face<D, M> s) const
{
  const NbrInfoBase* info = nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get();
  if (info == nullptr) {
    throw RuntimeError("PatchInfo::getFineNbrInfo: nbr_info is nullptr");
  }
  return dynamic_cast<const FineNbrInfo<M>&>(*info);
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
  requires is_valid_face<D, M>
inline bool
PatchInfo<D>::hasNbr(Face<D, M> s) const
{
  return nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()] != nullptr;
}

template<int D>
  requires is_supported_dimension<D>
inline bool
PatchInfo<D>::hasNbr() const
{
  return std::any_of(
    nbr_infos.begin(), nbr_infos.end(), [](const std::unique_ptr<NbrInfoBase>& nbr_info) {
      return nbr_info != nullptr;
    });
}

template<int D>
  requires is_supported_dimension<D> bool
PatchInfo<D>::hasCoarseParent() const
{
  return orth_on_parent != Orthant<D>::null();
}

template<int D>
  requires is_supported_dimension<D>
void
PatchInfo<D>::setNeighborLocalIndexes(const std::map<int, int>& id_to_local_index_map)
{
  for (size_t i = 0; i < nbr_infos.size(); i++) {
    if (nbr_infos[i] != nullptr) {
      nbr_infos[i]->setLocalIndexes(id_to_local_index_map);
    }
  }
}

template<int D>
  requires is_supported_dimension<D>
void
PatchInfo<D>::setNeighborGlobalIndexes(const std::map<int, int>& id_to_global_index_map)
{
  for (size_t i = 0; i < nbr_infos.size(); i++) {
    if (nbr_infos[i] != nullptr) {
      nbr_infos[i]->setGlobalIndexes(id_to_global_index_map);
    }
  }
}

template<int D>
  requires is_supported_dimension<D>
std::deque<int>
PatchInfo<D>::getNbrIds() const
{
  std::deque<int> retval;
  for (size_t i = 0; i < nbr_infos.size(); i++) {
    if (nbr_infos[i] != nullptr) {
      nbr_infos[i]->getNbrIds(retval);
    }
  }
  return retval;
}

template<int D>
  requires is_supported_dimension<D>
std::deque<int>
PatchInfo<D>::getNbrRanks() const
{
  std::deque<int> retval;
  for (size_t i = 0; i < nbr_infos.size(); i++) {
    if (nbr_infos[i] != nullptr) {
      nbr_infos[i]->getNbrRanks(retval);
    }
  }
  return retval;
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
void
PatchInfo<D>::serializeNeighbors(BufferWriter& writer) const
{
  std::bitset<Face<D, M>::number_of> has_nbr;
  for (Face<D, M> f : Face<D, M>::getValues()) {
    has_nbr[f.getIndex()] = hasNbr(f);
  }
  writer << has_nbr;
  for (Face<D, M> f : Face<D, M>::getValues()) {
    if (hasNbr(f)) {
      NbrType type = getNbrType(f);
      writer << type;
      switch (type) {
        case NbrType::Normal: {
          NormalNbrInfo<M> tmp = getNormalNbrInfo(f);
          writer << tmp;
        } break;
        case NbrType::Fine: {
          FineNbrInfo<M> tmp = getFineNbrInfo(f);
          writer << tmp;
        } break;
        case NbrType::Coarse: {
          CoarseNbrInfo<M> tmp = getCoarseNbrInfo(f);
          writer << tmp;
        } break;
        default:
          throw RuntimeError("Unsupported NbrType");
      }
    }
  }
  if constexpr (M > 0) {
    serializeNeighbors<M - 1>(writer);
  }
}

template<int D>
  requires is_supported_dimension<D>
template<int M>
void
PatchInfo<D>::deserializeNeighbors(BufferReader& reader)
{
  std::bitset<Face<D, M>::number_of> has_nbr;
  reader >> has_nbr;
  for (Face<D, M> f : Face<D, M>::getValues()) {
    if (has_nbr[f.getIndex()]) {
      NbrType type;
      reader >> type;
      NbrInfo<M>* info;
      switch (type) {
        case NbrType::Normal:
          info = new NormalNbrInfo<M>();
          reader >> *static_cast<NormalNbrInfo<M>*>(info);
          break;
        case NbrType::Fine:
          info = new FineNbrInfo<M>();
          reader >> *static_cast<FineNbrInfo<M>*>(info);
          break;
        case NbrType::Coarse:
          info = new CoarseNbrInfo<M>();
          reader >> *static_cast<CoarseNbrInfo<M>*>(info);
          break;
        default:
          throw RuntimeError("Unsupported NbrType");
      }

      setNbrInfo(f, info);
    }
  }
  if constexpr (M > 0) {
    deserializeNeighbors<M - 1>(reader);
  }
}

template<int D>
  requires is_supported_dimension<D>
int
PatchInfo<D>::serialize(char* buffer) const
{
  BufferWriter writer(buffer);
  writer << id;
  writer << ns;
  writer << refine_level;
  writer << parent_id;
  writer << parent_rank;
  writer << child_ids;
  writer << child_ranks;
  writer << rank;
  writer << orth_on_parent;
  writer << starts;
  writer << spacings;

  serializeNeighbors<D - 1>(writer);

  return writer.getPos();
}

template<int D>
  requires is_supported_dimension<D>
int
PatchInfo<D>::deserialize(char* buffer)
{
  BufferReader reader(buffer);
  reader >> id;
  reader >> ns;
  reader >> refine_level;
  reader >> parent_id;
  reader >> parent_rank;
  reader >> child_ids;
  reader >> child_ranks;
  reader >> rank;
  reader >> orth_on_parent;
  reader >> starts;
  reader >> spacings;

  deserializeNeighbors<D - 1>(reader);

  return reader.getPos();
}

// EXPLICIT INSTANTIATIONS

template class PatchInfo<2>;

template void
PatchInfo<2>::setNbrInfo(Face<2, 1> f, std::nullptr_t);
template void
PatchInfo<2>::setNbrInfo(Face<2, 0> f, std::nullptr_t);

template void
PatchInfo<2>::setNbrInfo(Face<2, 1> f, NbrInfo<1>* nbr_info);
template void
PatchInfo<2>::setNbrInfo(Face<2, 0> f, NbrInfo<0>* nbr_info);

template NbrType
PatchInfo<2>::getNbrType(Face<2, 1> f) const;
template NbrType
PatchInfo<2>::getNbrType(Face<2, 0> f) const;

template const NormalNbrInfo<1>&
PatchInfo<2>::getNormalNbrInfo(Face<2, 1> f) const;
template const NormalNbrInfo<0>&
PatchInfo<2>::getNormalNbrInfo(Face<2, 0> f) const;

template const CoarseNbrInfo<1>&
PatchInfo<2>::getCoarseNbrInfo(Face<2, 1> f) const;
template const CoarseNbrInfo<0>&
PatchInfo<2>::getCoarseNbrInfo(Face<2, 0> f) const;

template const FineNbrInfo<1>&
PatchInfo<2>::getFineNbrInfo(Face<2, 1> f) const;
template const FineNbrInfo<0>&
PatchInfo<2>::getFineNbrInfo(Face<2, 0> f) const;

template class PatchInfo<3>;

template void
PatchInfo<3>::setNbrInfo(Face<3, 2> f, std::nullptr_t);
template void
PatchInfo<3>::setNbrInfo(Face<3, 1> f, std::nullptr_t);
template void
PatchInfo<3>::setNbrInfo(Face<3, 0> f, std::nullptr_t);

template NbrType
PatchInfo<3>::getNbrType(Face<3, 2> f) const;
template NbrType
PatchInfo<3>::getNbrType(Face<3, 1> f) const;
template NbrType
PatchInfo<3>::getNbrType(Face<3, 0> f) const;

template const NormalNbrInfo<2>&
PatchInfo<3>::getNormalNbrInfo(Face<3, 2> f) const;
template const NormalNbrInfo<1>&
PatchInfo<3>::getNormalNbrInfo(Face<3, 1> f) const;
template const NormalNbrInfo<0>&
PatchInfo<3>::getNormalNbrInfo(Face<3, 0> f) const;

template const CoarseNbrInfo<2>&
PatchInfo<3>::getCoarseNbrInfo(Face<3, 2> f) const;
template const CoarseNbrInfo<1>&
PatchInfo<3>::getCoarseNbrInfo(Face<3, 1> f) const;
template const CoarseNbrInfo<0>&
PatchInfo<3>::getCoarseNbrInfo(Face<3, 0> f) const;

template const FineNbrInfo<2>&
PatchInfo<3>::getFineNbrInfo(Face<3, 2> f) const;
template const FineNbrInfo<1>&
PatchInfo<3>::getFineNbrInfo(Face<3, 1> f) const;
template const FineNbrInfo<0>&
PatchInfo<3>::getFineNbrInfo(Face<3, 0> f) const;

} // namespace ThunderEgg
