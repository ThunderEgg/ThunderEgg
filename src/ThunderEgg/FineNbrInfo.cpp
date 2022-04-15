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
#include <ThunderEgg/FineNbrInfo.h>

namespace ThunderEgg {

template<int D>
  requires is_supported_neighbor_dimension<D>
FineNbrInfo<D>::FineNbrInfo()
{
  ranks.fill(0);
  ids.fill(0);
  local_indexes.fill(-1);
  global_indexes.fill(-1);
}

template<int D>
  requires is_supported_neighbor_dimension<D>
FineNbrInfo<D>::~FineNbrInfo() = default;

template<int D>
  requires is_supported_neighbor_dimension<D>
FineNbrInfo<D>::FineNbrInfo(std::array<int, Orthant<D>::num_orthants> ids)
  : ids(ids)
{
  ranks.fill(0);
  local_indexes.fill(-1);
  global_indexes.fill(-1);
}

template<int D>
  requires is_supported_neighbor_dimension<D>
NbrType
FineNbrInfo<D>::getNbrType() const
{
  return NbrType::Fine;
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
FineNbrInfo<D>::getNbrIds(std::deque<int>& nbr_ids) const
{
  for (size_t i = 0; i < ids.size(); i++) {
    nbr_ids.push_back(ids[i]);
  }
};

template<int D>
  requires is_supported_neighbor_dimension<D>
void
FineNbrInfo<D>::getNbrRanks(std::deque<int>& nbr_ranks) const
{
  for (size_t i = 0; i < ranks.size(); i++) {
    nbr_ranks.push_back(ranks[i]);
  }
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
FineNbrInfo<D>::setGlobalIndexes(const std::map<int, int>& id_to_global_index_map)
{
  for (size_t i = 0; i < global_indexes.size(); i++) {
    global_indexes[i] = id_to_global_index_map.at(ids[i]);
  }
}

template<int D>
  requires is_supported_neighbor_dimension<D>
void
FineNbrInfo<D>::setLocalIndexes(const std::map<int, int>& id_to_local_index_map)
{
  for (size_t i = 0; i < local_indexes.size(); i++) {
    auto iter = id_to_local_index_map.find(ids[i]);
    if (iter != id_to_local_index_map.end()) {
      local_indexes[i] = iter->second;
    }
  }
}

template<int D>
  requires is_supported_neighbor_dimension<D>
int
FineNbrInfo<D>::serialize(char* buffer) const
{
  BufferWriter writer(buffer);
  writer << ranks;
  writer << ids;
  return writer.getPos();
}

template<int D>
  requires is_supported_neighbor_dimension<D>
int
FineNbrInfo<D>::deserialize(char* buffer)
{
  BufferReader reader(buffer);
  reader >> ranks;
  reader >> ids;
  return reader.getPos();
}

template<int D>
  requires is_supported_neighbor_dimension<D>
std::unique_ptr<NbrInfoBase>
FineNbrInfo<D>::clone() const
{
  return std::make_unique<FineNbrInfo<D>>(*this);
}

// EXPLICIT INSTANTIATIONS

template class FineNbrInfo<0>;
template class FineNbrInfo<1>;
template class FineNbrInfo<2>;

} // namespace ThunderEgg