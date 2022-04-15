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
#include <ThunderEgg/NormalNbrInfo.h>

namespace ThunderEgg {

template<int D>
  requires is_supported_neighbor_dimension<D>
NormalNbrInfo<D>::NormalNbrInfo()
{
}

template<int D>
  requires is_supported_neighbor_dimension<D>
NormalNbrInfo<D>::~NormalNbrInfo() = default;

template<int D>
  requires is_supported_neighbor_dimension<D>
NormalNbrInfo<D>::NormalNbrInfo(int id)
{
  this->id = id;
}

template<int D>
  requires is_supported_neighbor_dimension<D>
NbrType
NormalNbrInfo<D>::getNbrType() const
{
  return NbrType::Normal;
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
NormalNbrInfo<D>::getNbrIds(std::deque<int>& nbr_ids) const
{
  nbr_ids.push_back(id);
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
NormalNbrInfo<D>::getNbrRanks(std::deque<int>& nbr_ranks) const
{
  nbr_ranks.push_back(rank);
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
NormalNbrInfo<D>::setGlobalIndexes(const std::map<int, int>& id_to_global_index_map)
{
  global_index = id_to_global_index_map.at(id);
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
NormalNbrInfo<D>::setLocalIndexes(const std::map<int, int>& id_to_local_index_map)
{
  auto iter = id_to_local_index_map.find(id);
  if (iter != id_to_local_index_map.end()) {
    local_index = iter->second;
  }
}

template<int D>
  requires is_supported_neighbor_dimension<D>
int
NormalNbrInfo<D>::serialize(char* buffer) const
{
  BufferWriter writer(buffer);
  writer << rank;
  writer << id;
  return writer.getPos();
}

template<int D>
  requires is_supported_neighbor_dimension<D>
int
NormalNbrInfo<D>::deserialize(char* buffer)
{
  BufferReader reader(buffer);
  reader >> rank;
  reader >> id;
  return reader.getPos();
}

template<int D>
  requires is_supported_neighbor_dimension<D>
std::unique_ptr<NbrInfoBase>
NormalNbrInfo<D>::clone() const
{
  return std::make_unique<NormalNbrInfo<D>>(*this);
}

// EXPLICIT INSTANTIATIONS

template class NormalNbrInfo<0>;
template class NormalNbrInfo<1>;
template class NormalNbrInfo<2>;

} // namespace ThunderEgg