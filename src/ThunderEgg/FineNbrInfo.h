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
#ifndef THUNDEREGG_FINENBRINFO_H
#define THUNDEREGG_FINENBRINFO_H
/**
 * @file
 *
 * @brief FineNbrInfo class
 */
#include <ThunderEgg/BufferReader.h>
#include <ThunderEgg/BufferWriter.h>
#include <ThunderEgg/NbrInfo.h>
#include <ThunderEgg/Orthant.h>

namespace ThunderEgg {
/**
 * @brief Represents neighbors that are at a finer refinement level.
 *
 * @tparam D the number of Cartesian dimensions.
 */
template<int D>
  requires is_supported_neighbor_dimension<D>
class FineNbrInfo : public NbrInfo<D>
{
public:
  /**
   * @brief The mpi rank that the neighbor resides on.
   */
  std::array<int, Orthant<D>::num_orthants> ranks;
  /**
   * @brief The ids of the neighbors
   */
  std::array<int, Orthant<D>::num_orthants> ids;
  /**
   * @brief The global indexes of the neighbors
   */
  std::array<int, Orthant<D>::num_orthants> global_indexes;
  /**
   * @brief The local indexes of the neighbors
   */
  std::array<int, Orthant<D>::num_orthants> local_indexes;

  /**
   * @brief Construct a new empty FineNbrInfo object
   */
  FineNbrInfo();

  ~FineNbrInfo();

  /**
   * @brief Construct a new FineNbrInfo object
   *
   * @param ids the ids of the neighbors
   */
  FineNbrInfo(std::array<int, Orthant<D>::num_orthants> ids);

  NbrType
  getNbrType() const override;

  void
  getNbrIds(std::deque<int>& nbr_ids) const override;

  void
  getNbrRanks(std::deque<int>& nbr_ranks) const override;

  void
  setGlobalIndexes(const std::map<int, int>& id_to_global_index_map) override;

  void
  setLocalIndexes(const std::map<int, int>& id_to_local_index_map) override;

  int
  serialize(char* buffer) const override;

  int
  deserialize(char* buffer) override;

  std::unique_ptr<NbrInfoBase>
  clone() const override;
};

template<int D>
  requires is_supported_neighbor_dimension<D>
void
to_json(tpl::nlohmann::json& j, const FineNbrInfo<D>& n);

template<int D>
  requires is_supported_neighbor_dimension<D>
void
from_json(const tpl::nlohmann::json& j, FineNbrInfo<D>& n);

// EXPLICIT INSTANTIATIONS

extern template class FineNbrInfo<0>;
extern template class FineNbrInfo<1>;
extern template class FineNbrInfo<2>;

extern template void
to_json<0>(tpl::nlohmann::json& j, const FineNbrInfo<0>& n);
extern template void
to_json<1>(tpl::nlohmann::json& j, const FineNbrInfo<1>& n);
extern template void
to_json<2>(tpl::nlohmann::json& j, const FineNbrInfo<2>& n);

extern template void
from_json<0>(const tpl::nlohmann::json& j, FineNbrInfo<0>& n);
extern template void
from_json<1>(const tpl::nlohmann::json& j, FineNbrInfo<1>& n);
extern template void
from_json<2>(const tpl::nlohmann::json& j, FineNbrInfo<2>& n);

} // namespace ThunderEgg
#endif