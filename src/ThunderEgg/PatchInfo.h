/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2017-2021 Scott Aiton
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
#ifndef THUNDEREGG_PATCHINFO_H
#define THUNDEREGG_PATCHINFO_H
/**
 * @file
 *
 * @brief PatchInfo class
 */
#include <ThunderEgg/CoarseNbrInfo.h>
#include <ThunderEgg/Config.h>
#include <ThunderEgg/FineNbrInfo.h>
#include <ThunderEgg/NormalNbrInfo.h>
#include <ThunderEgg/Orthant.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Serializable.h>
#include <algorithm>
#include <bitset>

namespace ThunderEgg {
/**
 * @brief Contains metadata for a patch
 *
 * This contains metadata for a specific patch. Information like:
 * * The globally unique id of this patch
 * * The local and global indexes of this patch in the Domain
 * * The parent patch in the tree (if there is one)
 *
 * It also contains information for a patch's neighbor:
 * * What are the neighbors id?
 * * Are the neighbors at the same refinement level? Are coarser or finer?
 *
 * @tparam D the number of cartesian dimensions in the patch
 */
template<int D>
  requires is_supported_dimension<D>
class PatchInfo : public Serializable
{
private:
  /**
   * @brief Nbr info objects for each side.
   * If there is no neighbor, it should be set to nullptr.
   */
  std::array<std::unique_ptr<NbrInfoBase>, Face<D, D>::sum_of_faces> nbr_infos;

  /**
   * @brief Serialize neighbors into buffer
   *
   * @tparam M dimensionality of neighbor, this is recursive stops at 0
   * @param writer the buffer
   */
  template<int M>
  void
  serializeNeighbors(BufferWriter& writer) const;

  /**
   * @brief Deserialize neighbors from buffer
   *
   * @tparam M dimensionality of neighbor, this is recursive stops at 0
   * @param reader the buffer
   */
  template<int M>
  void
  deserializeNeighbors(BufferReader& reader);

public:
  /**
   * @brief The globally unique ID of the patch
   * This ID only needs to be unique within a Domain.
   */
  int id = 0;
  /**
   * @brief The local index of the patch in the Domain.
   */
  int local_index = 0;
  /**
   * @brief The global index of the patch in the Domain.
   */
  int global_index = 0;
  /**
   * @brief The refinement level
   */
  int refine_level = -1;
  /**
   * @brief The id of the parent patch.
   *
   * Set to -1 if there is no parent.
   */
  int parent_id = -1;
  /**
   * @brief the rank that the parent patch resides on
   */
  int parent_rank = -1;
  /**
   * @brief The id's of the children.
   *
   * Set to -1 if there are no children
   */
  std::array<int, Orthant<D>::num_orthants> child_ids;
  /**
   * @brief The ranks of the children
   *
   * Set to -1 if there are no children
   */
  std::array<int, Orthant<D>::num_orthants> child_ranks;
  /**
   * @brief Number of ghost cells on each side of the patch.
   */
  int num_ghost_cells = 0;
  /**
   * @brief MPI rank of this patch
   */
  int rank = -1;
  /**
   * @brief The orthant of the parent that this parent resides on.
   *
   * If the parent is the same size, it should be set to Orthant::null
   */
  Orthant<D> orth_on_parent = Orthant<D>::null();
  /**
   * @brief The number of cells in each direction
   */
  std::array<int, D> ns;
  /**
   * @brief The lower-left-bottom index of the patch
   */
  std::array<double, D> starts;
  /**
   * @brief The cell spacings in each direction
   */
  std::array<double, D> spacings;

  /**
   * @brief Construct a new Patch Info object
   * starts, ns, and spacings are all set to 0
   */
  PatchInfo();

  /**
   * @brief Copy constructor
   *
   * @param other_pinfo  object to copy
   */
  PatchInfo(const PatchInfo<D>& other_pinfo);

  /**
   * @brief Copy asisgnment
   *
   * @param other_pinfo the object to copy
   * @return PatchInfo<D>& this object
   */
  PatchInfo<D>&
  operator=(const PatchInfo<D>& other_pinfo);

  /**
   * @brief Compare the ids of the patches
   *
   * @param r right operand
   * @return true if r's id is lower
   * @return false if r's id is not lower
   */
  bool
  operator<(const PatchInfo<D>& r);

  /**
   * @brief Set the Nbr Info object to null
   *
   * @tparam M the dimensionality of the face
   * @param f the face
   */
  template<int M>
    requires is_valid_face<D, M>
  void
  setNbrInfo(Face<D, M> f, std::nullptr_t);

  /**
   * @brief Set the Nbr Info object on a given face
   *
   * @tparam M the dimensionality of the face
   * @param s the face
   * @param nbr_info the neighbor info object, this patchinfo will take ownership of it
   */
  template<int M>
    requires is_valid_face<D, M>
  void
  setNbrInfo(Face<D, M> f, NbrInfo<M>* nbr_info);

  /**
   * @brief Get the NbrType for a side
   *
   * @param s the side
   * @return The NbrType
   */
  template<int M>
    requires is_valid_face<D, M>
  NbrType
  getNbrType(Face<D, M> s) const;

  /**
   * @brief Get the NormalNbrInfo object for a side
   *
   * Neighbor must be of Normal type, otherwise behavior is undefined.
   *
   * @param s the side
   * @return const NormalNbrInfo<D>& the object
   */
  template<int M>
    requires is_valid_face<D, M>
  const NormalNbrInfo<M>&
  getNormalNbrInfo(Face<D, M> s) const;

  /**
   * @brief Get the CoarseNbrInfo object
   *
   * @param s the side
   * @return const CoarseNbrInfo<D>& the object
   */
  template<int M>
    requires is_valid_face<D, M>
  const CoarseNbrInfo<M>&
  getCoarseNbrInfo(Face<D, M> s) const;

  /**
   * @brief Get the FineNbrInfo object
   *
   * Neighbor must be of Fine type, otherwise behavior is undefined.
   *
   * @param s the side
   * @return const FineNbrInfo<D>& the object
   */
  template<int M>
    requires is_valid_face<D, M>
  const FineNbrInfo<M>&
  getFineNbrInfo(Face<D, M> s) const;

  /**
   * @brief Return whether the patch has a neighbor
   *
   * @param s the side
   * @return true if the is neighbor
   * @return false if at domain boundary
   */
  template<int M>
    requires is_valid_face<D, M> bool
  hasNbr(Face<D, M> s) const;

  /**
   * @brief Return if this patch has a neighbor
   */
  bool
  hasNbr() const;

  /**
   * @brief Return whether the patch has a coarser parent
   *
   * @return true if there is a parent
   */
  bool
  hasCoarseParent() const;

  /**
   * @brief Set the local indexes in the NbrInfo objects
   *
   * @param id_to_local_index_map map from id to local_index
   */
  void
  setNeighborLocalIndexes(const std::map<int, int>& id_to_local_index_map);

  /**
   * @brief Set the global indexes in the NbrInfo objects
   *
   * @param id_to_global_index_map map form id to global_index
   */
  void
  setNeighborGlobalIndexes(const std::map<int, int>& id_to_global_index_map);

  /**
   * @brief return a vector of neighbor ids
   */
  std::deque<int>
  getNbrIds() const;

  /**
   * @brief return a vector of neighbor ranks
   */
  std::deque<int>
  getNbrRanks() const;

  /**
   * @brief Serialize the object from a buffer
   *
   * @param buffer the buffer (can be nullptr for dry run to get the size)
   * @return int the number of bytes written
   */
  int
  serialize(char* buffer) const;

  /**
   * @brief Deserialize the object from a buffer
   *
   * @param buffer the buffer
   * @return int the number of bytes read
   */
  int
  deserialize(char* buffer);
};

template<int D>
  requires is_supported_dimension<D>
void
to_json(tpl::nlohmann::json& j, const PatchInfo<D>& pinfo);

template<int D>
  requires is_supported_dimension<D>
void
from_json(const tpl::nlohmann::json& j, PatchInfo<D>& pinfo);

extern template class PatchInfo<2>;
extern template class PatchInfo<3>;

extern template void
to_json(tpl::nlohmann::json& j, const PatchInfo<2>& pinfo);
extern template void
to_json(tpl::nlohmann::json& j, const PatchInfo<3>& pinfo);

extern template void
from_json(const tpl::nlohmann::json& j, PatchInfo<2>& pinfo);
extern template void
from_json(const tpl::nlohmann::json& j, PatchInfo<3>& pinfo);

} // namespace ThunderEgg
#endif