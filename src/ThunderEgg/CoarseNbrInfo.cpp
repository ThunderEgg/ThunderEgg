/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#include <ThunderEgg/CoarseNbrInfo.h>

namespace ThunderEgg {

template<int D>
  requires is_supported_neighbor_dimension<D>
CoarseNbrInfo<D>::CoarseNbrInfo() = default;
template<int D>
  requires is_supported_neighbor_dimension<D>
CoarseNbrInfo<D>::~CoarseNbrInfo() = default;

template<int D>
  requires is_supported_neighbor_dimension<D>
CoarseNbrInfo<D>::CoarseNbrInfo(int id, Orthant<D> orth_on_coarse)
{
  this->id = id;
  this->orth_on_coarse = orth_on_coarse;
}

template<int D>
  requires is_supported_neighbor_dimension<D>
NbrType
CoarseNbrInfo<D>::getNbrType() const
{
  return NbrType::Coarse;
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
CoarseNbrInfo<D>::getNbrIds(std::deque<int>& nbr_ids) const
{
  nbr_ids.push_back(id);
};

template<int D>
  requires is_supported_neighbor_dimension<D>
void
CoarseNbrInfo<D>::getNbrRanks(std::deque<int>& nbr_ranks) const
{
  nbr_ranks.push_back(rank);
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
CoarseNbrInfo<D>::setGlobalIndexes(const std::map<int, int>& id_to_global_index_map)
{
  global_index = id_to_global_index_map.at(id);
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
CoarseNbrInfo<D>::setLocalIndexes(const std::map<int, int>& id_to_local_index_map)
{
  auto iter = id_to_local_index_map.find(id);
  if (iter != id_to_local_index_map.end()) {
    local_index = iter->second;
  }
}

template<int D>
  requires is_supported_neighbor_dimension<D>
int
CoarseNbrInfo<D>::serialize(char* buffer) const
{
  BufferWriter writer(buffer);
  writer << rank;
  writer << id;
  writer << orth_on_coarse;
  return writer.getPos();
}

template<int D>
  requires is_supported_neighbor_dimension<D>
int
CoarseNbrInfo<D>::deserialize(char* buffer)
{
  BufferReader reader(buffer);
  reader >> rank;
  reader >> id;
  reader >> orth_on_coarse;
  return reader.getPos();
}

template<int D>
  requires is_supported_neighbor_dimension<D>
std::unique_ptr<NbrInfoBase>
CoarseNbrInfo<D>::clone() const
{
  return std::make_unique<CoarseNbrInfo<D>>(*this);
}

// EXPLICIT INSTANTIATIONS

template class CoarseNbrInfo<0>;
template class CoarseNbrInfo<1>;
template class CoarseNbrInfo<2>;

} // namespace ThunderEgg