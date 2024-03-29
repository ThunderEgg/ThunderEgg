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

#ifndef THUNDEREGG_SCHUR_InterfaceDomain_H
#define THUNDEREGG_SCHUR_InterfaceDomain_H
/**
 * @file
 *
 * @brief InterfaceDomain class
 */
#include <ThunderEgg/Schur/Interface.h>
#include <ThunderEgg/Schur/PatchIfaceInfo.h>
#include <ThunderEgg/Vector.h>
#include <deque>
namespace ThunderEgg {
namespace Schur {
/**
 * @brief Represents the Schur compliment domain of the problem.
 *
 * This class mainly manages a set of interfaces that makes up the Schur compliment system. It is
 * responsible for setting up the indexing of the interfaces, which is used in the rest of the
 * ThunderEgg library.
 *
 * This class will set a global index for each interface, and will provide three sets of local
 * indexes.
 *
 * The patch_iface set of indexes are the local indexes for a vector of interfaces that line up with
 * each piinfo object.
 *
 * The row set of indexes are for a row distributed parallel matrix (where a rank contains the
 * entire row of the matrix)
 *
 * The column set of indexes are for a column distributed parallel matrix.
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
class InterfaceDomain
{
private:
  std::array<int, D - 1> iface_ns;
  Domain<D> domain;

  /**
   * @brief Vector of PatchIfaceInfo pointers where index in the vector corresponds to the patch's
   * local index
   */
  std::vector<std::shared_ptr<const PatchIfaceInfo<D>>> piinfos;
  /**
   * @brief Vector of Interfaces pointers where index in the vector corresponds to the
   * interfaces's local index
   */
  std::vector<std::shared_ptr<const Interface<D>>> interfaces;
  /**
   * @brief The global number of interfaces
   */
  int num_global_ifaces = 0;

  /**
   * @brief Index all of column, row, and patch interface local indexes for the interface system
   *
   * @param id_to_iface_map map of Interface id to Interface objects
   * @param piinfos the vector PatchIfaceInfo objects for this processor
   * @param interfaces (output) this will be updated with a new vector of Interface objects, the
   * position in the vector cooresponds to the Interface's local index
   */
  static void IndexIfacesLocal(const std::map<int, std::shared_ptr<Interface<D>>>& id_to_iface_map,
                               const std::vector<std::shared_ptr<PatchIfaceInfo<D>>>& piinfos,
                               std::vector<std::shared_ptr<Interface<D>>>& interfaces)
  {
    int curr_local_index = 0;

    // index interface objects first
    interfaces.clear();
    interfaces.reserve(id_to_iface_map.size());
    for (auto pair : id_to_iface_map) {
      auto iface = pair.second;

      iface->local_index = curr_local_index;
      interfaces.push_back(iface);

      for (auto patch : iface->patches) {
        if (patch.type.isNormal() || patch.type.isFineToFine() || patch.type.isCoarseToCoarse()) {
          patch.getNonConstPiinfo()->getIfaceInfo(patch.side)->patch_local_index = curr_local_index;
          patch.getNonConstPiinfo()->getIfaceInfo(patch.side)->col_local_index = curr_local_index;
          patch.getNonConstPiinfo()->getIfaceInfo(patch.side)->row_local_index = curr_local_index;
        } else if (patch.type.isFineToCoarse()) {
          patch.getNonConstPiinfo()->getCoarseIfaceInfo(patch.side)->coarse_col_local_index =
            curr_local_index;
        } else if (patch.type.isCoarseToFine()) {
          auto iface_info = patch.getNonConstPiinfo()->getFineIfaceInfo(patch.side);
          for (size_t i = 0; i < iface_info->fine_col_local_indexes.size(); i++) {
            if (iface_info->fine_ids[i] == iface->id) {
              iface_info->fine_col_local_indexes[i] = curr_local_index;
              break;
            }
          }
        }
      }

      curr_local_index++;
    }

    // perform rest of necessary indexing
    IndexRemainingColIfacesLocal(curr_local_index, interfaces);
    IndexRemainingRowIfacesLocal(curr_local_index, interfaces);
    IndexRemainingPatchIfacesLocal(curr_local_index, piinfos);
  }
  /**
   * @brief Get the remaining unset column local indexes in the given PatchIfaceInfo object
   *
   * @param piinfo the PatchIfaceInfo object
   * @param id_to_local_indexes_to_set map from id to pointer to local index to set
   */
  static void GetRemainginColIfacesLocalForPatch(
    std::shared_ptr<PatchIfaceInfo<D>> piinfo,
    std::map<int, std::set<int*>>& id_to_local_indexes_to_set)
  {
    for (Side<D> s : Side<D>::getValues()) {
      if (piinfo->pinfo.hasNbr(s)) {
        auto iface_info = piinfo->getIfaceInfo(s);
        if (iface_info->col_local_index == -1) {
          id_to_local_indexes_to_set[iface_info->id].insert(&iface_info->col_local_index);
        }

        NbrType nbr_type = piinfo->pinfo.getNbrType(s);

        if (nbr_type == NbrType::Coarse) {
          auto coarse_iface_info = piinfo->getCoarseIfaceInfo(s);
          if (coarse_iface_info->coarse_col_local_index == -1) {
            id_to_local_indexes_to_set[coarse_iface_info->coarse_id].insert(
              &coarse_iface_info->coarse_col_local_index);
          }
        } else if (nbr_type == NbrType::Fine) {
          auto fine_iface_info = piinfo->getFineIfaceInfo(s);
          for (size_t i = 0; i < fine_iface_info->fine_col_local_indexes.size(); i++) {
            if (fine_iface_info->fine_col_local_indexes[i] == -1) {
              id_to_local_indexes_to_set[fine_iface_info->fine_ids[i]].insert(
                &fine_iface_info->fine_col_local_indexes[i]);
            }
          }
        }
      }
    }
  }
  /**
   * @brief Index the remaining unset column local indexes
   *
   * @param curr_local_index the current index
   * @param interfaces the vector Interface objects
   */
  static void IndexRemainingColIfacesLocal(
    int curr_local_index,
    const std::vector<std::shared_ptr<Interface<D>>>& interfaces)
  {
    std::map<int, std::set<int*>> id_to_local_indexes_to_set;
    for (auto iface : interfaces) {
      for (auto patch : iface->patches) {
        auto piinfo = patch.getNonConstPiinfo();

        if (patch.type.isNormal() || patch.type.isFineToFine() || patch.type.isCoarseToCoarse()) {
          GetRemainginColIfacesLocalForPatch(piinfo, id_to_local_indexes_to_set);
        }
      }
    }
    for (auto id_and_local_indexes : id_to_local_indexes_to_set) {
      for (int* local_index_to_set : id_and_local_indexes.second) {
        *local_index_to_set = curr_local_index;
      }
      curr_local_index++;
    }
  }
  /**
   * @brief Index the remaining unset row local indexes
   *
   * @param curr_local_index the current index
   * @param interfaces the vector Interface objects
   */
  static void IndexRemainingRowIfacesLocal(
    int curr_local_index,
    const std::vector<std::shared_ptr<Interface<D>>>& interfaces)
  {
    std::map<int, std::set<int*>> id_to_local_indexes_to_set;
    for (auto iface : interfaces) {
      for (auto patch : iface->patches) {
        auto piinfo = patch.getNonConstPiinfo();

        for (Side<D> s : Side<D>::getValues()) {
          if (piinfo->pinfo.hasNbr(s)) {
            auto iface_info = piinfo->getIfaceInfo(s);

            if (iface_info->row_local_index == -1) {
              id_to_local_indexes_to_set[iface_info->id].insert(&iface_info->row_local_index);
            }
          }
        }
      }
    }
    for (auto id_and_local_indexes : id_to_local_indexes_to_set) {
      for (int* local_index_to_set : id_and_local_indexes.second) {
        *local_index_to_set = curr_local_index;
      }
      curr_local_index++;
    }
  }
  /**
   * @brief Index the remaining unset patch interface local indexes
   *
   * @param curr_local_index the current index
   * @param interfaces the vector PatchIfaceInfo objects
   */
  static void IndexRemainingPatchIfacesLocal(
    int curr_local_index,
    const std::vector<std::shared_ptr<PatchIfaceInfo<D>>>& piinfos)
  {
    for (auto piinfo : piinfos) {
      for (Side<D> s : Side<D>::getValues()) {
        if (piinfo->pinfo.hasNbr(s)) {
          auto iface_info = piinfo->getIfaceInfo(s);

          if (iface_info->patch_local_index == -1) {
            iface_info->patch_local_index = curr_local_index;
            curr_local_index++;
          }
        }
      }
    }
  }
  /**
   * @brief Set global indexes for all of the interfaces, local indexes should already be set
   *
   * @param interfaces the vector of Interface objects for this processor
   * @param piinfos the vector PatchIfaceInfo objects for this processor
   */
  static void IndexIfacesGlobal(const std::vector<std::shared_ptr<Interface<D>>>& interfaces,
                                const std::vector<std::shared_ptr<PatchIfaceInfo<D>>>& piinfos)
  {
    // get starting global index for this rank
    int starting_global_index;
    int num_local_interfaces = (int)interfaces.size();
    MPI_Scan(&num_local_interfaces, &starting_global_index, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    starting_global_index -= num_local_interfaces;

    // index local interfaces first
    for (auto iface : interfaces) {
      iface->global_index = starting_global_index + iface->local_index;

      for (auto patch : iface->patches) {
        if (patch.type.isNormal() || patch.type.isFineToFine() || patch.type.isCoarseToCoarse()) {
          auto iface_info = patch.getNonConstPiinfo()->getIfaceInfo(patch.side);
          iface_info->global_index = iface->global_index;
        } else if (patch.type.isFineToCoarse()) {
          auto iface_info = patch.getNonConstPiinfo()->getCoarseIfaceInfo(patch.side);
          iface_info->coarse_global_index = iface->global_index;
        } else if (patch.type.isCoarseToFine()) {
          auto iface_info = patch.getNonConstPiinfo()->getFineIfaceInfo(patch.side);
          for (size_t i = 0; i < iface_info->fine_col_local_indexes.size(); i++) {
            if (iface_info->fine_ids[i] == iface->id) {
              iface_info->fine_global_indexes[i] = iface->global_index;
              break;
            }
          }
        }
      }
    }
    SendAndReceiveGlobalIndexes(interfaces, piinfos);
  }
  /**
   * @brief Do the necessary communication to get the global indexes from other processors
   *
   * @param interfaces the set of Interface objects
   * @param piinfos  the set of PatchIfaceInfo objects
   */
  static void SendAndReceiveGlobalIndexes(
    const std::vector<std::shared_ptr<Interface<D>>>& interfaces,
    const std::vector<std::shared_ptr<PatchIfaceInfo<D>>>& piinfos)
  {
    std::map<int, std::map<int, std::set<int*>>> rank_to_id_to_global_indexes_to_set;
    std::deque<std::vector<int>> recv_buffers;
    std::vector<MPI_Request> recv_requests;
    SetupGlobalIndexRecvRequests(
      interfaces, piinfos, rank_to_id_to_global_indexes_to_set, recv_buffers, recv_requests);

    std::deque<std::vector<int>> send_buffers;
    std::vector<MPI_Request> send_requests;
    SetupGlobalIndexSendRequests(interfaces, send_buffers, send_requests);

    size_t num_recvs = recv_requests.size();
    for (size_t i = 0; i < num_recvs; i++) {
      int index;
      MPI_Status status;
      MPI_Waitany(recv_requests.size(), recv_requests.data(), &index, &status);
      const std::vector<int>& buffer = recv_buffers[index];
      const auto& id_to_global_indexes_to_set =
        rank_to_id_to_global_indexes_to_set[status.MPI_SOURCE];

      // set the global indexes
      size_t curr_index = 0;
      for (auto pair : id_to_global_indexes_to_set) {
        for (int* global_index_to_set : pair.second) {
          *global_index_to_set = buffer[curr_index];
        }
        curr_index++;
      }
    }

    MPI_Waitall((int)send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
  }
  /**
   * @brief Setup the MPI_Irecv calls
   *
   * @param interfaces the vector of Interface objects
   * @param piinfos the vector of PatchIfaceInfo objects
   * @param rank_to_id_to_global_indexes_to_set (output) Map from rank of incoming process to
   * id of interface to pointers to global index values that have to be set
   * @param recv_buffers (output) the buffers for the recv requests. When recv is finished,
   * global index values will by sorted by the interface's id.
   * @param recv_requests (output) MPI_Irecv request status, one for each incoming rank.
   */
  static void SetupGlobalIndexRecvRequests(
    const std::vector<std::shared_ptr<Interface<D>>>& interfaces,
    const std::vector<std::shared_ptr<PatchIfaceInfo<D>>>& piinfos,
    std::map<int, std::map<int, std::set<int*>>>& rank_to_id_to_global_indexes_to_set,
    std::deque<std::vector<int>>& recv_buffers,
    std::vector<MPI_Request>& recv_requests)
  {
    GetGlobalIndexesToSet(interfaces, piinfos, rank_to_id_to_global_indexes_to_set);

    for (auto pair : rank_to_id_to_global_indexes_to_set) {
      recv_buffers.emplace_back(pair.second.size());
      MPI_Request request;
      MPI_Irecv(recv_buffers.back().data(),
                (int)recv_buffers.back().size(),
                MPI_INT,
                pair.first,
                0,
                MPI_COMM_WORLD,
                &request);
      recv_requests.push_back(request);
    }
  }
  /**
   * @brief Get pointers to the global indexes that have to be set for a given patch
   *
   * @param piinfo the PatchIfaceInfo object
   * @param rank_to_id_to_global_indexes_to_set (output) Map from rank of incoming process to
   * id of interface to pointers to global index values that have to be set
   */
  static void GetGlobalIndexesToSetForOuterInterfaces(
    std::shared_ptr<PatchIfaceInfo<D>> piinfo,
    std::map<int, std::map<int, std::set<int*>>>& rank_to_id_to_global_indexes_to_set)
  {
    for (Side<D> s : Side<D>::getValues()) {
      if (piinfo->pinfo.hasNbr(s)) {
        NbrType nbr_type = piinfo->pinfo.getNbrType(s);

        if (nbr_type == NbrType::Coarse) {
          auto coarse_iface_info = piinfo->getCoarseIfaceInfo(s);
          if (coarse_iface_info->coarse_global_index == -1) {
            rank_to_id_to_global_indexes_to_set[coarse_iface_info->coarse_rank]
                                               [coarse_iface_info->coarse_id]
                                                 .insert(&coarse_iface_info->coarse_global_index);
          }
        } else if (nbr_type == NbrType::Fine) {
          auto fine_iface_info = piinfo->getFineIfaceInfo(s);
          for (size_t i = 0; i < fine_iface_info->fine_global_indexes.size(); i++) {
            if (fine_iface_info->fine_global_indexes[i] == -1) {
              rank_to_id_to_global_indexes_to_set[fine_iface_info->fine_ranks[i]]
                                                 [fine_iface_info->fine_ids[i]]
                                                   .insert(
                                                     &fine_iface_info->fine_global_indexes[i]);
            }
          }
        }
      }
    }
  }
  /**
   * @brief Get pointers to the global indexes that have to be set.
   *
   * @param interfaces the vector of Interface objects
   * @param piinfos the vector of PatchIfaceInfo objects
   * @param rank_to_id_to_global_indexes_to_set (output) Map from rank of incoming process to
   * id of interface to pointers to global index values that have to be set
   */
  static void GetGlobalIndexesToSet(
    const std::vector<std::shared_ptr<Interface<D>>>& interfaces,
    const std::vector<std::shared_ptr<PatchIfaceInfo<D>>>& piinfos,
    std::map<int, std::map<int, std::set<int*>>>& rank_to_id_to_global_indexes_to_set)
  {
    // patch interfaces
    GetGlobalIndexesToSetForPatchInterfaces(piinfos, rank_to_id_to_global_indexes_to_set);

    // rows and columns
    for (auto iface : interfaces) {
      for (auto patch : iface->patches) {
        auto piinfo = patch.getNonConstPiinfo();

        // row
        for (Side<D> s : Side<D>::getValues()) {
          // the inner interfaces of this patch will affect this interface
          if (piinfo->pinfo.hasNbr(s)) {
            auto iface_info = piinfo->getIfaceInfo(s);

            if (iface_info->global_index == -1) {
              rank_to_id_to_global_indexes_to_set[iface_info->rank][iface_info->id].insert(
                &iface_info->global_index);
            }
          }
        }
        // column
        if (patch.type.isNormal() || patch.type.isFineToFine() || patch.type.isCoarseToCoarse()) {
          // this interface will affect the values of the outer interfaces
          GetGlobalIndexesToSetForOuterInterfaces(piinfo, rank_to_id_to_global_indexes_to_set);
        }
      }
    }
  }
  /**
   * @brief Get the global indexes to set For PatchInterfacesInfo objects
   *
   * @param piinfos the vector of PatchIfaceInfo objects
   * @param rank_to_id_to_global_indexes_to_set (output) Map from rank of incoming process to
   * id of interface to pointers to global index values that have to be set
   */
  static void GetGlobalIndexesToSetForPatchInterfaces(
    const std::vector<std::shared_ptr<PatchIfaceInfo<D>>>& piinfos,
    std::map<int, std::map<int, std::set<int*>>>& rank_to_id_to_global_indexes_to_set)
  {
    for (auto piinfo : piinfos) {
      for (Side<D> s : Side<D>::getValues()) {
        if (piinfo->pinfo.hasNbr(s)) {
          auto iface_info = piinfo->getIfaceInfo(s);

          if (iface_info->global_index == -1) {
            rank_to_id_to_global_indexes_to_set[iface_info->rank][iface_info->id].insert(
              &iface_info->global_index);
          }
        }
      }
    }
  }
  /**
   * @brief Setup the MPI_Isend calls
   *
   * @param interfaces the vector of Interface objects
   * @param send_buffers (output) the buffers for the send requests. Should not be deallocated
   * until sends are done.
   * @param send_requests (output) MPI_Isend request status, one for each outgoing rank
   */
  static void SetupGlobalIndexSendRequests(
    const std::vector<std::shared_ptr<Interface<D>>>& interfaces,
    std::deque<std::vector<int>>& send_buffers,
    std::vector<MPI_Request>& send_requests)
  {
    std::map<int, std::set<std::pair<int, int>>> rank_to_id_and_global_index_pairs;
    GetGlobalIndexesToSend(interfaces, rank_to_id_and_global_index_pairs);

    for (auto pair : rank_to_id_and_global_index_pairs) {
      send_buffers.emplace_back();
      std::vector<int>& send_buffer = send_buffers.back();
      send_buffer.reserve(pair.second.size());

      for (auto id_global_index_pair : pair.second) {
        send_buffer.push_back(id_global_index_pair.second);
      }

      MPI_Request request;
      MPI_Isend(send_buffer.data(),
                (int)send_buffer.size(),
                MPI_INT,
                pair.first,
                0,
                MPI_COMM_WORLD,
                &request);
      send_requests.push_back(request);
    }
  }
  /**
   * @brief Get the global indexes that have to be sent
   *
   * @param interfaces the vector of Interface objects
   * @param rank_to_id_and_global_index_pairs (output) Map from rank of incoming process to set
   * of pairs of ids and global indexes
   */
  static void GetGlobalIndexesToSend(
    const std::vector<std::shared_ptr<Interface<D>>>& interfaces,
    std::map<int, std::set<std::pair<int, int>>>& rank_to_id_and_global_index_pairs)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (auto iface : interfaces) {
      for (auto patch : iface->patches) {
        auto piinfo = patch.piinfo;
        // patch interfaces
        if (piinfo->pinfo.rank != rank) {
          rank_to_id_and_global_index_pairs[piinfo->pinfo.rank].emplace(iface->id,
                                                                        iface->global_index);
        }
        // cols
        for (Side<D> s : Side<D>::getValues()) {
          // the inner interfaces of this patch will affect this interface
          if (piinfo->pinfo.hasNbr(s)) {
            auto iface_info = piinfo->getIfaceInfo(s);

            if (iface_info->rank != rank) {
              rank_to_id_and_global_index_pairs[iface_info->rank].emplace(iface->id,
                                                                          iface->global_index);
            }
          }
        }
        // row
        if (patch.type.isNormal() || patch.type.isCoarseToCoarse() || patch.type.isFineToFine()) {
          // this interface will affect the values of the outer interfaces
          GetGlobalIndexesToSendForOuterInterfaces(
            iface, piinfo, rank_to_id_and_global_index_pairs);
        }
      }
    }
  }
  /**
   * @brief Get global indexes that have to be sent for patch
   *
   * @param interface the Interface that we are sending the global index from
   * @param piinfo the PatchIfaceInfo object
   * @param rank_to_id_to_global_indexes_to_set (output) Map from rank of incoming process to
   * id of interface to pointers to global index values that have to be set
   */
  static void GetGlobalIndexesToSendForOuterInterfaces(
    std::shared_ptr<const Interface<D>> interface,
    std::shared_ptr<const PatchIfaceInfo<D>> piinfo,
    std::map<int, std::set<std::pair<int, int>>>& rank_to_id_and_global_index_pairs)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (Side<D> s : Side<D>::getValues()) {
      if (piinfo->pinfo.hasNbr(s)) {
        NbrType nbr_type = piinfo->pinfo.getNbrType(s);

        if (nbr_type == NbrType::Coarse) {
          auto coarse_iface_info = piinfo->getCoarseIfaceInfo(s);
          if (coarse_iface_info->coarse_rank != rank) {
            rank_to_id_and_global_index_pairs[coarse_iface_info->coarse_rank].emplace(
              interface->id, interface->global_index);
          }
        } else if (nbr_type == NbrType::Fine) {
          auto fine_iface_info = piinfo->getFineIfaceInfo(s);
          for (size_t i = 0; i < fine_iface_info->fine_col_local_indexes.size(); i++) {
            if (fine_iface_info->fine_ranks[i] != rank) {
              rank_to_id_and_global_index_pairs[fine_iface_info->fine_ranks[i]].emplace(
                interface->id, interface->global_index);
            }
          }
        }
      }
    }
  }

public:
  InterfaceDomain() = default;
  /**
   * @brief Create a InterfaceDomain from a given Domain
   *
   * @param domain the Domain
   */
  explicit InterfaceDomain(const Domain<D>& domain)
    : domain(domain)
  {
    iface_ns.fill(domain.getNs()[0]);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::shared_ptr<PatchIfaceInfo<D>>> piinfos_non_const;
    piinfos.reserve(domain.getNumLocalPatches());
    piinfos_non_const.reserve(domain.getNumLocalPatches());
    for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
      piinfos_non_const.emplace_back(new PatchIfaceInfo<D>(pinfo));
      piinfos.push_back(piinfos_non_const.back());
    }

    std::map<int, std::map<int, std::shared_ptr<Schur::Interface<D>>>> rank_id_iface_map;
    std::vector<std::shared_ptr<Schur::PatchIfaceInfo<D>>> off_proc_piinfos;
    Interface<D>::EnumerateIfacesFromPiinfoVector(piinfos, rank_id_iface_map, off_proc_piinfos);

    std::vector<std::shared_ptr<Schur::Interface<D>>> interfaces_non_const;
    IndexIfacesLocal(rank_id_iface_map[rank], piinfos_non_const, interfaces_non_const);
    IndexIfacesGlobal(interfaces_non_const, piinfos_non_const);

    interfaces.reserve(interfaces_non_const.size());
    for (auto iface : interfaces_non_const) {
      interfaces.push_back(iface);
    }

    int num_ifaces = interfaces.size();
    MPI_Allreduce(&num_ifaces, &num_global_ifaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  /**
   * @brief Get the number of local Interfaces on this rank
   *
   * @return int the number of local Interfaces
   */
  int getNumLocalInterfaces() const { return interfaces.size(); }
  /**
   * @brief Get the number of Interfaces on all ranks
   *
   * @return int the number of Interfaces on all ranks
   */
  int getNumGlobalInterfaces() const { return num_global_ifaces; }
  /**
   * @brief Get the vector Interfaces objects for this rank
   *
   * The location of each Interface in the vector will coorespond to the Interafce's local
   * index
   *
   * @return const std::vector<std::shared_ptr<const Interface<D>>> the vector of
   * Interface objects
   */
  const std::vector<std::shared_ptr<const Interface<D>>> getInterfaces() const
  {
    return interfaces;
  }
  /**
   * @brief Get the vector PatchIfaceInfo objects for this rank
   *
   * The location of each PatchIfaceInfo in the vector will coorespond to the patch's
   * local index
   *
   * @return const std::vector<std::shared_ptr<const PatchIfaceInfo<D>>>& the vector of
   * PatchIfaceInfo objects
   */
  const std::vector<std::shared_ptr<const PatchIfaceInfo<D>>>& getPatchIfaceInfos() const
  {
    return piinfos;
  }
  /**
   * @brief Get the Domain object that cooresponds to this InterfaceDomain
   *
   * @return std::shared_ptr<Domain<D>> the Domain object
   */
  const Domain<D>& getDomain() const { return domain; }

  /**
   * @brief Get a new vector for the schur compliment system
   *
   * @return Vector<D - 1> the vector
   */
  Vector<D - 1> getNewVector() const
  {
    return Vector<D - 1>(domain.getCommunicator(), iface_ns, 1, getNumLocalInterfaces(), 0);
  }
};
extern template class InterfaceDomain<2>;
extern template class InterfaceDomain<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
