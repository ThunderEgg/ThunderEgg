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
class PatchInfo : public Serializable
{
private:
  /**
   * @brief Nbr info objects for each side.
   * If there is no neighbor, it should be set to nullptr.
   */
  std::array<std::unique_ptr<NbrInfoBase>, Face<D, D>::sum_of_faces> nbr_infos;

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
  PatchInfo()
  {
    starts.fill(0);
    ns.fill(1);
    spacings.fill(1);
    child_ids.fill(-1);
    child_ranks.fill(-1);
  }

  /**
   * @brief Copy constructor
   *
   * @param other_pinfo  object to copy
   */
  PatchInfo(const PatchInfo<D>& other_pinfo)
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
  /**
   * @brief Copy asisgnment
   *
   * @param other_pinfo the object to copy
   * @return PatchInfo<D>& this object
   */
  PatchInfo<D>& operator=(const PatchInfo<D>& other_pinfo)
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
  /**
   * @brief Compare the ids of the patches
   *
   * @param l left operand
   * @param r right operand
   * @return true if r's id is lower
   * @return false if r's id is not lower
   */
  friend bool operator<(const PatchInfo& l, const PatchInfo& r) { return l.id < r.id; }
  template<int M>
  void setNbrInfo(Face<D, M> f, std::nullptr_t)
  {
    nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()] = nullptr;
  }
  /**
   * @brief Set the Nbr Info object on a given face
   *
   * @tparam M the dimensionality of the face
   * @param s the face
   * @param nbr_info the neighbor info object, this patchinfo will take ownership of it
   */
  template<int M>
  void setNbrInfo(Face<D, M> f, NbrInfo<M>* nbr_info)
  {
    nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()].reset(nbr_info);
  }
  /**
   * @brief Get the NbrType for a side
   *
   * @param s the side
   * @return The NbrType
   */
  template<int M>
  NbrType getNbrType(Face<D, M> s) const
  {
    return nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()]->getNbrType();
  }
  /**
   * @brief Get the NormalNbrInfo object for a side
   *
   * Neighbor must be of Normal type, otherwise behavior is undefined.
   *
   * @param s the side
   * @return NormalNbrInfo<D>& the object
   */
  template<int M>
  NormalNbrInfo<M>& getNormalNbrInfo(Face<D, M> s) const
  {
    return *dynamic_cast<NormalNbrInfo<M>*>(
      nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get());
  }
  /**
   * @brief Get the CoarseNbrInfo object
   *

   * @param s the side
   * @return CoarseNbrInfo<D>& the object
  */
  template<int M>
  CoarseNbrInfo<M>& getCoarseNbrInfo(Face<D, M> s) const
  {
    return *dynamic_cast<CoarseNbrInfo<M>*>(
      nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get());
  }
  /**
   * @brief Get the FineNbrInfo object
   *
   * Neighbor must be of Fine type, otherwise behavior is undefined.
   *
   * @param s the side
   * @return FineNbrInfo<D>& the object
   */
  template<int M>
  FineNbrInfo<M>& getFineNbrInfo(Face<D, M> s) const
  {
    return *dynamic_cast<FineNbrInfo<M>*>(nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get());
  }
  /**
   * @brief Return whether the patch has a neighbor
   *
   * @param s the side
   * @return true if the is neighbor
   * @return false if at domain boundary
   */
  template<int M>
  inline bool hasNbr(Face<D, M> s) const
  {
    return nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()] != nullptr;
  }
  /**
   * @brief Return if this patch has a neighbor
   */
  inline bool hasNbr() const
  {
    return std::any_of(
      nbr_infos.begin(), nbr_infos.end(), [](const std::unique_ptr<NbrInfoBase>& nbr_info) {
        return nbr_info != nullptr;
      });
  }
  /**
   * @brief Return whether the patch has a coarser parent
   *
   * @return true if there is a parent
   */
  inline bool hasCoarseParent() const { return orth_on_parent != Orthant<D>::null(); }
  /**
   * @brief Set the local indexes in the NbrInfo objects
   *
   * @param id_to_local_index_map map from id to local_index
   */
  void setNeighborLocalIndexes(const std::map<int, int>& id_to_local_index_map)
  {
    for (size_t i = 0; i < nbr_infos.size(); i++) {
      if (nbr_infos[i] != nullptr) {
        nbr_infos[i]->setLocalIndexes(id_to_local_index_map);
      }
    }
  }
  /**
   * @brief Set the global indexes in the NbrInfo objects
   *
   * @param id_to_global_index_map map form id to global_index
   */
  void setNeighborGlobalIndexes(const std::map<int, int>& id_to_global_index_map)
  {
    for (size_t i = 0; i < nbr_infos.size(); i++) {
      if (nbr_infos[i] != nullptr) {
        nbr_infos[i]->setGlobalIndexes(id_to_global_index_map);
      }
    }
  }
  /**
   * @brief return a vector of neighbor ids
   */
  std::deque<int> getNbrIds() const
  {
    std::deque<int> retval;
    for (size_t i = 0; i < nbr_infos.size(); i++) {
      if (nbr_infos[i] != nullptr) {
        nbr_infos[i]->getNbrIds(retval);
      }
    }
    return retval;
  }
  /**
   * @brief return a vector of neighbor ranks
   */
  std::deque<int> getNbrRanks() const
  {
    std::deque<int> retval;
    for (size_t i = 0; i < nbr_infos.size(); i++) {
      if (nbr_infos[i] != nullptr) {
        nbr_infos[i]->getNbrRanks(retval);
      }
    }
    return retval;
  }
  template<int M>
  void serializeNeighbors(BufferWriter& writer) const
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
  template<int M>
  void deserializeNeighbors(BufferReader& reader)
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
  int serialize(char* buffer) const
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
  int deserialize(char* buffer)
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
};

template<int D>
void
to_json(tpl::nlohmann::json& j, const PatchInfo<D>& pinfo);

template<int D>
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