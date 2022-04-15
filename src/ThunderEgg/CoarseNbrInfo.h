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
#ifndef THUNDEREGG_COARSENBRINFO_H
#define THUNDEREGG_COARSENBRINFO_H
/**
 * @file
 *
 * @brief CoarseNbrInfo class
 */
#include <ThunderEgg/BufferReader.h>
#include <ThunderEgg/BufferWriter.h>
#include <ThunderEgg/NbrInfo.h>
#include <ThunderEgg/Orthant.h>

namespace ThunderEgg {
/**
 * @brief Represents a neighbor that is at a coarser refinement level.
 *
 * @tparam D the number of Cartesian dimensions.
 */
template<int D>
  requires is_supported_neighbor_dimension<D>
class CoarseNbrInfo : public NbrInfo<D>
{
public:
  /**
   * @brief The mpi rank that the neighbor resides on.
   */
  int rank = 0;
  /**
   * @brief The id of the neighbor
   */
  int id = 0;
  /**
   * @brief The local index of the neighbor
   */
  int local_index = -1;
  /**
   * @brief The global index of the neighbor
   */
  int global_index = -1;
  /**
   * @brief The orthant that this patch in relation to the coarser patch's interface.
   */
  Orthant<D> orth_on_coarse;
  /**
   * @brief Construct a new empty CoarseNbrInfo object
   */
  CoarseNbrInfo();

  ~CoarseNbrInfo();
  /**
   * @brief Construct a new CoarseNbrInfo object
   *
   * @param id the id of the neighbor
   * @param orth_on_coarse The orthant of the neighboring patch's interface that the this patch
   * lies along.
   */
  CoarseNbrInfo(int id, Orthant<D> orth_on_coarse);

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
to_json(tpl::nlohmann::json& j, const CoarseNbrInfo<D>& n);

template<int D>
  requires is_supported_neighbor_dimension<D>
void
from_json(const tpl::nlohmann::json& j, CoarseNbrInfo<D>& n);

// EXPLICIT INSTANTIATIONS

extern template class CoarseNbrInfo<0>;
extern template class CoarseNbrInfo<1>;
extern template class CoarseNbrInfo<2>;

extern template void
to_json<0>(tpl::nlohmann::json& j, const CoarseNbrInfo<0>& n);
extern template void
to_json<1>(tpl::nlohmann::json& j, const CoarseNbrInfo<1>& n);
extern template void
to_json<2>(tpl::nlohmann::json& j, const CoarseNbrInfo<2>& n);

extern template void
from_json<0>(const tpl::nlohmann::json& j, CoarseNbrInfo<0>& n);
extern template void
from_json<1>(const tpl::nlohmann::json& j, CoarseNbrInfo<1>& n);
extern template void
from_json<2>(const tpl::nlohmann::json& j, CoarseNbrInfo<2>& n);

} // namespace ThunderEgg
#endif