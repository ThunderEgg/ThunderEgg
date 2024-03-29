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

#ifndef THUNDEREGG_MPIGHOSTFILLER_H
#define THUNDEREGG_MPIGHOSTFILLER_H
/**
 * @file
 *
 * @brief MPIGhostFiller class
 */

#include <ThunderEgg/DimensionalArray.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/GhostFillingType.h>

namespace ThunderEgg {
/**
 * @brief Parallell ghostfiller implimented with MPI
 *
 * There are three functions that have to be overridden in derived classes.
 * fillGhostCellsForNbrPatch, fillGhostCellsForEdgeNbrPatch, fillGhostCellsForCornerNbrPatch, and
 * fillGhostCellsForLocalPatch
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
class MPIGhostFiller : public GhostFiller<D>
{
private:
  /**
   * @brief Structure fore remote call information used during construction.
   *
   * @tparam M the dimension of the face that the call is on
   */
  template<int M>
  class RemoteCallPrototype
  {
  public:
    /**
     * @brief The id of the neighboring patch
     */
    int id;
    /**
     * @brief The face of the local patch that the neighboring patch is on
     */
    Face<D, M> face;
    /**
     * @brief The type of neighbor
     */
    NbrType nbr_type;
    /**
     * @brief The orthant on the on the coarse patch, or null
     */
    Orthant<M> orthant;
    /**
     * @brief The local index of the local patch the we are ghost filling from
     */
    int local_index;
    /**
     * @brief Construct a new Remote Call Prototype object
     *
     * @param id the id of the patch
     * @param face the face of the neighboring patch is on
     * @param nbr_type the type of neighbor
     * @param orthant the orthant on the coarse patch, or null
     * @param local_index the lcoal index of the patch that is being filled from
     */
    RemoteCallPrototype(int id,
                        Face<D, M> face,
                        NbrType nbr_type,
                        Orthant<M> orthant,
                        int local_index)
      : id(id)
      , face(face)
      , nbr_type(nbr_type)
      , orthant(orthant)
      , local_index(local_index)
    {}
    /**
     * @brief Sort by id, face.opposite(), nbr_type, orthant, and then local_index
     *
     * @param other the other object to compare to
     * @return the comparison
     */
    bool operator<(const RemoteCallPrototype& other) const
    {
      return std::forward_as_tuple(id, face.opposite(), nbr_type, orthant, local_index) <
             std::forward_as_tuple(
               other.id, other.face.opposite(), other.nbr_type, other.orthant, other.local_index);
    }
  };

  /**
   * @brief Incoming ghost structure used in construction
   *
   * @tparam M the dimension of the face
   */
  template<int M>
  class IncomingGhostPrototype
  {
  public:
    /**
     * @brief The id of the patch that is being filled
     */
    int id;
    /**
     * @brief The face of the patch that ghosts are on
     */
    Face<D, M> face;
    /**
     * @brief the local index of the patch being filled
     */
    int local_index;
    /**
     * @brief Construct a new Incoming Ghost Prototype object
     *
     * @param id the id of the patch
     * @param face the face of the patch that the ghosts are on
     * @param local_index the local index of the patch
     */
    IncomingGhostPrototype(int id, Face<D, M> face, int local_index)
      : id(id)
      , face(face)
      , local_index(local_index)
    {}
    /**
     * @brief Sort by id, face, and then local_index
     *
     * @param other the other object to compare to
     * @return the comparison
     */
    bool operator<(const IncomingGhostPrototype& other) const
    {
      return std::forward_as_tuple(id, face, local_index) <
             std::forward_as_tuple(other.id, other.face, other.local_index);
    }
  };

  /**
   * @brief Struct for remote call
   *
   * @tparam M the dimension of the face
   */
  template<int M>
  struct RemoteCall
  {
    /**
     * @brief The side of the patch that the neighboring patch is on
     */
    Face<D, M> face;
    /**
     * @brief The type of neighbor
     */
    NbrType nbr_type;
    /**
     * @brief The orthant on the coarse patch, or null
     */
    Orthant<M> orthant;
    /**
     * @brief The local index of the patch that we are filling from
     *
     */
    int local_index;
    /**
     * @brief offset in buffer outgoing buffer
     */
    size_t offset;
    /**
     * @brief Construct a new Remote Call object
     *
     * @param prototype the protype call
     * @param offset the offset in the buffer
     */
    RemoteCall(const RemoteCallPrototype<M>& prototype, size_t offset)
      : face(prototype.face)
      , nbr_type(prototype.nbr_type)
      , orthant(prototype.orthant)
      , local_index(prototype.local_index)
      , offset(offset)
    {}
  };

  /**
   * @brief The struct for incoming ghosts
   *
   * @tparam M the dimension of the face
   */
  template<int M>
  struct IncomingGhost
  {
    /**
     * @brief The side of the patch that the ghosts are on
     */
    Face<D, M> face;
    /**
     * @brief The local index of the patch that is being filled
     */
    int local_index;
    /**
     * @brief The offset in the incoming buffer
     */
    size_t offset;
    /**
     * @brief Construct a new Incoming Ghost object
     *
     * @param prototype the prototype object
     * @param offset the offset in the incoming buffer
     */
    IncomingGhost(const IncomingGhostPrototype<M>& prototype, size_t offset)
      : face(prototype.face)
      , local_index(prototype.local_index)
      , offset(offset)
    {}
  };

  /**
   * @brief Struct for local call
   *
   * @tparam M the dimension of the face
   */
  template<int M>
  struct LocalCall
  {
    /**
     * @brief The face that the neighbor is on
     */
    Face<D, M> face;
    /**
     * @brief The neighbor type
     */
    NbrType nbr_type;
    /**
     * @brief The orthant on the coarse patch
     */
    Orthant<M> orthant;
    /**
     * @brief The local index of the patch that is being filled from
     */
    int local_index;
    /**
     * @brief The local index of the patch that is being filled
     */
    int nbr_local_index;
    /**
     * @brief Construct a new Local Call object
     *
     * @param face the face that the loca call is on
     * @param nbr_type the neighbor type
     * @param orthant the orthant on the coarse patch, or null
     * @param local_index the local index of the patch being filled from
     * @param nbr_local_index the local index of the patch being filled
     */
    LocalCall(Face<D, M> face,
              NbrType nbr_type,
              Orthant<M> orthant,
              int local_index,
              size_t nbr_local_index)
      : face(face)
      , nbr_type(nbr_type)
      , orthant(orthant)
      , local_index(local_index)
      , nbr_local_index(nbr_local_index)
    {}
  };

  /**
   * @brief Data for a remote call
   */
  struct RemoteCallSet
  {
    /**
     * @brief The rank that we are communicating with
     */
    int rank;
    /**
     * @brief the length of the send buffer
     */
    size_t send_buffer_length = 0;
    /**
     * @brief deque of RemoteCall with dimension M
     *
     * @tparam M the dimension
     */
    template<int M>
    using RemoteCallDeque = std::deque<RemoteCall<M>>;
    /**
     * @brief An array of RemoteCall deques
     */
    DimensionalArray<D, RemoteCallDeque> remote_calls;
    /**
     * @brief the length of the recv buffer
     */
    size_t recv_buffer_length = 0;
    /**
     * @brief deque of IncomingGhost with dimension M
     *
     * @tparam M the dimension
     */
    template<int M>
    using IncomingGhostDeque = std::deque<IncomingGhost<M>>;
    /**
     * @brief An array of IncomingGhost deques
     */
    DimensionalArray<D, IncomingGhostDeque> incoming_ghosts;

    /**
     * @brief Construct a new Remote Call Set object
     *
     * @param rank the rank that is being communicated with
     */
    explicit RemoteCallSet(int rank)
      : rank(rank)
    {}
  };

  /**
   * @brief The vector of RemoteCallSet. one for each rank being communicated with
   */
  std::vector<RemoteCallSet> remote_call_sets;

  /**
   * @brief deque of LocalCall with dimension M
   *
   * @tparam M the dimension
   */
  template<int M>
  using LocalCallDeque = std::deque<LocalCall<M>>;
  /**
   * @brief array of deques of local calls to be made
   */
  DimensionalArray<D, LocalCallDeque> local_calls;

  /**
   * @brief Information need for representing ghosts in the buffer as full View object.
   *
   * @tparam M the dimension of the face
   */
  template<int M>
  class GhostViewInfo
  {
  private:
    /**
     * @brief The strides for each face value
     */
    std::array<std::array<int, D + 1>, Face<D, M>::number_of> strides;
    /**
     * @brief The ghost start for each face value
     */
    std::array<std::array<int, D + 1>, Face<D, M>::number_of> ghost_start;
    /**
     * @brief The start for each face value
     */
    std::array<std::array<int, D + 1>, Face<D, M>::number_of> start;
    /**
     * @brief The end for each face value
     */
    std::array<std::array<int, D + 1>, Face<D, M>::number_of> end;
    /**
     * @brief The ghost end for each face value
     */
    std::array<std::array<int, D + 1>, Face<D, M>::number_of> ghost_end;
    /**
     * @brief The sizes for each face value
     */
    std::array<size_t, Face<D, M>::number_of> sizes;

  public:
    /**
     * @brief Construct a new Ghost Local Data Info object
     */
    GhostViewInfo() = default;
    /**
     * @brief Construct a new Ghost Local Data Info object for a given Domain
     *
     * @param domain the domain
     */
    explicit GhostViewInfo(const Domain<D>& domain)
    {
      std::array<int, D + 1> ns;
      for (int i = 0; i < D; i++) {
        ns[i] = domain.getNs()[i];
      }
      ns[D] = 1;
      int num_ghost_cells = domain.getNumGhostCells();
      for (Face<D, M> face : Face<D, M>::getValues()) {
        std::array<int, D + 1> ghost_ns = ns;

        std::array<int, D + 1> face_ghost_start;
        face_ghost_start.fill(0);

        std::array<int, D + 1> face_start;
        face_start.fill(0);

        std::array<int, D + 1> face_end;
        face_end = ns;
        for (int& v : face_end) {
          v--;
        }

        std::array<int, D + 1> face_ghost_end;
        face_ghost_end = ns;
        for (int& v : face_ghost_end) {
          v--;
        }

        for (Side<D> side : face.getSides()) {
          ghost_ns[side.getAxisIndex()] = num_ghost_cells;
          if (side.isLowerOnAxis()) {
            face_ghost_start[side.getAxisIndex()] = -num_ghost_cells;
            face_ghost_end[side.getAxisIndex()] = -1;
            face_end[side.getAxisIndex()] = -1;
          } else {
            face_ghost_start[side.getAxisIndex()] = ns[side.getAxisIndex()];
            face_start[side.getAxisIndex()] = ns[side.getAxisIndex()];
            face_ghost_end[side.getAxisIndex()] = ns[side.getAxisIndex()] + num_ghost_cells - 1;
          }
        }
        strides[face.getIndex()][0] = 1;
        for (size_t i = 1; i < D + 1; i++) {
          strides[face.getIndex()][i] = ghost_ns[i - 1] * strides[face.getIndex()][i - 1];
        }
        sizes[face.getIndex()] = ghost_ns[D - 1] * strides[face.getIndex()][D - 1];
        int offset = 0;
        for (Side<D> side : face.getSides()) {
          if (side.isHigherOnAxis()) {
            offset += -(ns[side.getAxisIndex()] + num_ghost_cells) *
                      strides[face.getIndex()][side.getAxisIndex()];
          }
        }
        ghost_start[face.getIndex()] = face_ghost_start;
        start[face.getIndex()] = face_start;
        end[face.getIndex()] = face_end;
        ghost_end[face.getIndex()] = face_ghost_end;
      }
    }

    /**
     * @brief Get the the buffer size for ghosts on a given face
     *
     * @param face the face
     * @return size_t the size
     */
    size_t getSize(Face<D, M> face) const { return sizes[face.getIndex()]; }

    PatchView<const double, D> getPatchView(double* buffer_ptr,
                                            Face<D, M> face,
                                            int num_components) const
    {
      std::array<int, D + 1> my_end = end[face.getIndex()];
      std::array<int, D + 1> my_ghost_end = ghost_end[face.getIndex()];
      my_end[D] = num_components - 1;
      my_ghost_end[D] = num_components - 1;
      return PatchView<const double, D>(buffer_ptr,
                                        strides[face.getIndex()],
                                        ghost_start[face.getIndex()],
                                        start[face.getIndex()],
                                        my_end,
                                        my_ghost_end);
    }
  };

  /**
   * @brief GhostLocalDatInfo for each face dimension
   */
  DimensionalArray<D, GhostViewInfo> ghost_local_data_infos;

  /**
   * @brief The domain that this ghostfiller operates on
   */
  Domain<D> domain;

  /**
   * @brief the fill type to use
   */
  GhostFillingType fill_type;

  /**
   * @brief Get the View object for the buffer
   *
   * @param buffer_ptr pointer to the ghost cells position in the buffer
   * @param side  the side that the ghost cells are on
   * @return PatchView<D> the View object
   */
  template<int M>
  PatchView<const double, D> getPatchViewForBuffer(double* buffer_ptr,
                                                   Face<D, M> face,
                                                   int num_components) const
  {
    const GhostViewInfo<M>& gld_info = ghost_local_data_infos.template get<M>();
    return gld_info.getPatchView(buffer_ptr, face, num_components);
  }
  /**
   * @brief Post send requests
   *
   * @param buffers the allocated buffers
   * @return std::vector<MPI_Request> the send requests
   */
  std::vector<MPI_Request> postRecvs(std::vector<std::vector<double>>& buffers) const
  {
    std::vector<MPI_Request> recv_requests(buffers.size());
    for (size_t i = 0; i < recv_requests.size(); i++) {
      MPI_Irecv(buffers[i].data(),
                buffers[i].size(),
                MPI_DOUBLE,
                remote_call_sets[i].rank,
                0,
                MPI_COMM_WORLD,
                &recv_requests[i]);
    }
    return recv_requests;
  }

  /**
   * @brief A the values from the buffer to the ghost cells of patches
   *
   * @tparam M the dimension of the face
   * @param remote_call_set the remote calls
   * @param buffer the buffer
   * @param u the vector that the ghost values are being filled for
   */
  template<int M>
  void addRecvBufferToGhost(const RemoteCallSet& remote_call_set,
                            std::vector<double>& buffer,
                            const Vector<D>& u) const
  {
    for (const IncomingGhost<M>& incoming_ghost :
         remote_call_set.incoming_ghosts.template get<M>()) {
      int local_index = incoming_ghost.local_index;
      Face<D, M> face = incoming_ghost.face;
      size_t buffer_offset = incoming_ghost.offset;

      PatchView<const double, D> local_view = u.getPatchView(local_index);
      double* buffer_ptr = buffer.data() + buffer_offset * u.getNumComponents();
      PatchView<const double, D> buffer_view =
        getPatchViewForBuffer(buffer_ptr, face, u.getNumComponents());
      std::array<size_t, D - M> start;
      start.fill(0);
      std::array<size_t, D - M> end;
      end.fill(domain.getNumGhostCells() - 1);
      Loop::Nested<D - M>(start, end, [&](const std::array<size_t, D - M>& offset) {
        View<double, M + 1> local_slice = local_view.getGhostSliceOn(face, offset);
        View<double, M + 1> buffer_slice = buffer_view.getGhostSliceOn(face, offset);
        Loop::OverInteriorIndexes<M + 1>(local_slice, [&](const std::array<int, M + 1>& coord) {
          local_slice[coord] += buffer_slice[coord];
        });
      });
    }
  }
  /**
   * @brief process recv requests when they are ready
   *
   * @param requests the recv requests
   * @param buffers the recv buffers
   * @param u the vector to fill ghost values in
   */
  void processRecvs(std::vector<MPI_Request>& requests,
                    std::vector<std::vector<double>>& buffers,
                    const Vector<D>& u) const
  {
    size_t num_requests = requests.size();
    for (size_t i = 0; i < num_requests; i++) {
      int finished_index;
      MPI_Waitany(requests.size(), requests.data(), &finished_index, MPI_STATUS_IGNORE);
      switch (fill_type) {
        case GhostFillingType::Corners:
          if constexpr (D >= 2) {
            addRecvBufferToGhost<0>(remote_call_sets[finished_index], buffers[finished_index], u);
          }
        case GhostFillingType::Edges:
          if constexpr (D == 3) {
            addRecvBufferToGhost<1>(remote_call_sets[finished_index], buffers[finished_index], u);
          }
        case GhostFillingType::Faces:
          addRecvBufferToGhost<D - 1>(remote_call_sets[finished_index], buffers[finished_index], u);
          break;
        default:
          throw RuntimeError("Unsupported GhostFilling Type");
      }
    }
  }

  /**
   * @brief Fill the ghost values in the send buffer
   *
   * @tparam M the dimension of the face
   * @param remote_call_set the remote calls to be made
   * @param buffer the buffer that is being fileld
   * @param u the vector that the ghost values are being filled for
   */
  template<int M>
  void fillSendBuffer(const RemoteCallSet& remote_call_set,
                      std::vector<double>& buffer,
                      const Vector<D>& u) const
  {
    for (const RemoteCall<M>& call : remote_call_set.remote_calls.template get<M>()) {
      const PatchInfo<D>& pinfo = domain.getPatchInfoVector()[call.local_index];
      Face<D, M> face = call.face.opposite();
      PatchView<const double, D> local_view = u.getPatchView(call.local_index);
      double* buffer_ptr = buffer.data() + call.offset * u.getNumComponents();

      // create View objects for the buffer
      PatchView<const double, D> buffer_view =
        getPatchViewForBuffer(buffer_ptr, face, u.getNumComponents());

      // make the call
      fillGhostCellsForNbrPatchPriv(
        pinfo, local_view, buffer_view, call.face, call.nbr_type, call.orthant);
    }
  }

  /**
   * @brief fill buffers and post send requests
   *
   * @param buffers the allocated buffers
   * @param u the vector to fill buffers from
   * @return std::vector<MPI_Request> the requests
   */
  std::vector<MPI_Request> postSends(std::vector<std::vector<double>>& buffers,
                                     const Vector<D>& u) const
  {
    std::vector<MPI_Request> send_requests(remote_call_sets.size());
    for (size_t i = 0; i < remote_call_sets.size(); i++) {
      switch (fill_type) {
        case GhostFillingType::Corners:
          if constexpr (D >= 2) {
            fillSendBuffer<0>(remote_call_sets[i], buffers[i], u);
          }
        case GhostFillingType::Edges:
          if constexpr (D == 3) {
            fillSendBuffer<1>(remote_call_sets[i], buffers[i], u);
          }
        case GhostFillingType::Faces:
          fillSendBuffer<D - 1>(remote_call_sets[i], buffers[i], u);
          break;
        default:
          throw RuntimeError("Unsupported GhostFilling Type");
      }
      MPI_Isend(buffers[i].data(),
                buffers[i].size(),
                MPI_DOUBLE,
                remote_call_sets[i].rank,
                0,
                MPI_COMM_WORLD,
                &send_requests[i]);
    }
    return send_requests;
  }

  /**
   * @brief Fill the ghost cells for the neighboring patch
   *
   * @param pinfo the patch that ghost cells are being filled from
   * @param local_view the view for patch that ghost cells are being filled from
   * @param nbr_view  the view for the neighboring patch, where ghost cells are being
   * filled.
   * @param face the face that the neighboring patch is on
   * @param nbr_type the type of neighbor
   * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is
   * coarser
   */
  template<int M>
  void fillGhostCellsForNbrPatchPriv(const PatchInfo<D>& pinfo,
                                     const PatchView<const double, D>& local_view,
                                     const PatchView<const double, D>& nbr_view,
                                     Face<D, M> face,
                                     NbrType nbr_type,
                                     Orthant<M> orthant_on_coarse) const
  {
    if constexpr (M == D - 1) {
      fillGhostCellsForNbrPatch(pinfo, local_view, nbr_view, face, nbr_type, orthant_on_coarse);
    } else if constexpr (M > 0) {
      fillGhostCellsForEdgeNbrPatch(pinfo, local_view, nbr_view, face, nbr_type, orthant_on_coarse);
    } else {
      fillGhostCellsForCornerNbrPatch(pinfo, local_view, nbr_view, face, nbr_type);
    }
  }

  /**
   * @brief Process local fill calls for a given face dimension.
   *
   * The will make the call for dimension M-1
   *
   * @tparam M the dimension
   * @param u the vector that the ghost are being filled for
   */
  template<int M>
  void processLocalFills(const Vector<D>& u) const
  {
    for (const LocalCall<M>& call : local_calls.template get<M>()) {
      const PatchInfo<D>& pinfo = domain.getPatchInfoVector()[call.local_index];
      PatchView<const double, D> local_view = u.getPatchView(call.local_index);
      PatchView<const double, D> nbr_view = u.getPatchView(call.nbr_local_index);
      fillGhostCellsForNbrPatchPriv(
        pinfo, local_view, nbr_view, call.face, call.nbr_type, call.orthant);
    }
    if constexpr (M > 0) {
      processLocalFills<M - 1>(u);
    }
  }

  /**
   * @brief Add a normal neighbor call for given dimension
   *
   * @tparam M the dimension of the face
   * @param pinfo the PatchInfo object for the patch that is being filled from
   * @param f the the face
   * @param my_local_calls the local calls deque
   * @param rank_to_remote_call_prototypes the remote calls map
   * @param rank_to_incoming_ghost_prototypes the incoming ghosts map
   */
  template<int M>
  static void addNormalNbrCalls(
    const PatchInfo<D>& pinfo,
    Face<D, M> f,
    std::deque<LocalCall<M>>& my_local_calls,
    std::map<int, std::set<RemoteCallPrototype<M>>>& rank_to_remote_call_prototypes,
    std::map<int, std::set<IncomingGhostPrototype<M>>>& rank_to_incoming_ghost_prototypes)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto nbrinfo = pinfo.getNormalNbrInfo(f);
    if (nbrinfo.rank == rank) {
      my_local_calls.emplace_back(
        f, NbrType::Normal, Orthant<M>::null(), pinfo.local_index, nbrinfo.local_index);
    } else {
      rank_to_remote_call_prototypes[nbrinfo.rank].emplace(
        nbrinfo.id, f, NbrType::Normal, Orthant<M>::null(), pinfo.local_index);
      rank_to_incoming_ghost_prototypes[nbrinfo.rank].emplace(pinfo.id, f, pinfo.local_index);
    }
  }

  /**
   * @brief Add a fine neighbor call for given dimension
   *
   * @tparam M the dimension of the face
   * @param pinfo the PatchInfo object for the patch that is being filled from
   * @param f the the face
   * @param my_local_calls the local calls deque
   * @param rank_to_remote_call_prototypes the remote calls map
   * @param rank_to_incoming_ghost_prototypes the incoming ghosts map
   */
  template<int M>
  static void addFineNbrCalls(
    const PatchInfo<D>& pinfo,
    Face<D, M> f,
    std::deque<LocalCall<M>>& my_local_calls,
    std::map<int, std::set<RemoteCallPrototype<M>>>& rank_to_remote_call_prototypes,
    std::map<int, std::set<IncomingGhostPrototype<M>>>& rank_to_incoming_ghost_prototypes)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto nbrinfo = pinfo.getFineNbrInfo(f);
    for (size_t i = 0; i < Orthant<M>::num_orthants; i++) {
      if (nbrinfo.ranks[i] == rank) {
        my_local_calls.emplace_back(
          f, NbrType::Fine, Orthant<M>(i), pinfo.local_index, nbrinfo.local_indexes[i]);
      } else {
        rank_to_remote_call_prototypes[nbrinfo.ranks[i]].emplace(
          nbrinfo.ids[i], f, NbrType::Fine, Orthant<M>(i), pinfo.local_index);
        rank_to_incoming_ghost_prototypes[nbrinfo.ranks[i]].emplace(pinfo.id, f, pinfo.local_index);
      }
    }
  }

  /**
   * @brief Add a coarse neighbor call for given dimension
   *
   * @tparam M the dimension of the face
   * @param pinfo the PatchInfo object for the patch that is being filled from
   * @param f the the face
   * @param my_local_calls the local calls deque
   * @param rank_to_remote_call_prototypes the remote calls map
   * @param rank_to_incoming_ghost_prototypes the incoming ghosts map
   */
  template<int M>
  static void addCoarseNbrCalls(
    const PatchInfo<D>& pinfo,
    Face<D, M> f,
    std::deque<LocalCall<M>>& my_local_calls,
    std::map<int, std::set<RemoteCallPrototype<M>>>& rank_to_remote_call_prototypes,
    std::map<int, std::set<IncomingGhostPrototype<M>>>& rank_to_incoming_ghost_prototypes)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto nbrinfo = pinfo.getCoarseNbrInfo(f);
    if (nbrinfo.rank == rank) {
      my_local_calls.emplace_back(
        f, NbrType::Coarse, nbrinfo.orth_on_coarse, pinfo.local_index, nbrinfo.local_index);
    } else {
      rank_to_remote_call_prototypes[nbrinfo.rank].emplace(
        nbrinfo.id, f, NbrType::Coarse, nbrinfo.orth_on_coarse, pinfo.local_index);
      rank_to_incoming_ghost_prototypes[nbrinfo.rank].emplace(pinfo.id, f, pinfo.local_index);
    }
  }

  /**
   * @brief Enumerate calls on a given face dimension
   *
   * @tparam M the dimension of the faces
   * @param my_local_calls the local calls
   * @param rank_to_remote_call_sets a map from rank to RemoteCallSet objects
   */
  template<int M>
  void enumerateCalls(std::deque<LocalCall<M>>& my_local_calls,
                      std::map<int, RemoteCallSet>& rank_to_remote_call_sets) const
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::map<int, std::set<RemoteCallPrototype<M>>> rank_to_remote_call_prototypes;
    std::map<int, std::set<IncomingGhostPrototype<M>>> rank_to_incoming_ghost_prototypes;

    for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
      for (Face<D, M> f : Face<D, M>::getValues()) {
        if (pinfo.hasNbr(f)) {
          switch (pinfo.getNbrType(f)) {
            case NbrType::Normal:
              addNormalNbrCalls(pinfo,
                                f,
                                my_local_calls,
                                rank_to_remote_call_prototypes,
                                rank_to_incoming_ghost_prototypes);
              break;
            case NbrType::Fine:
              addFineNbrCalls(pinfo,
                              f,
                              my_local_calls,
                              rank_to_remote_call_prototypes,
                              rank_to_incoming_ghost_prototypes);
              break;
            case NbrType::Coarse:
              addCoarseNbrCalls(pinfo,
                                f,
                                my_local_calls,
                                rank_to_remote_call_prototypes,
                                rank_to_incoming_ghost_prototypes);
              break;
            default:
              throw RuntimeError("Unsupported NbrType");
          }
        }
      }
    }
    const GhostViewInfo<M>& ghost_local_data_info = ghost_local_data_infos.template get<M>();
    for (const auto& pair : rank_to_remote_call_prototypes) {
      int rank = pair.first;
      rank_to_remote_call_sets.emplace(rank, rank);
      RemoteCallSet& remote_call_set = rank_to_remote_call_sets.at(rank);
      std::tuple<int, Face<D, M>> prev_id_side;
      for (const RemoteCallPrototype<M>& call : pair.second) {
        size_t offset = remote_call_set.send_buffer_length;

        // calculate length in buffer need for ghost cells
        size_t length = ghost_local_data_info.getSize(call.face);
        // add length to buffer length
        // if its the same side of the patch, resuse the previous buffer space
        if (std::make_tuple(call.id, call.face) == prev_id_side) {
          offset -= length;
        } else {
          remote_call_set.send_buffer_length += length;
        }
        remote_call_set.remote_calls.template get<M>().emplace_back(call, offset);
        prev_id_side = std::make_tuple(call.id, call.face);
      }
    }
    for (const auto& pair : rank_to_incoming_ghost_prototypes) {
      RemoteCallSet& remote_call_set = rank_to_remote_call_sets.at(pair.first);
      for (const IncomingGhostPrototype<M>& prototype : pair.second) {
        // add length for ghosts to buffer length
        size_t length = ghost_local_data_info.getSize(prototype.face);
        size_t offset = remote_call_set.recv_buffer_length;
        remote_call_set.recv_buffer_length += length;

        // add ghost to incoming ghosts
        remote_call_set.incoming_ghosts.template get<M>().emplace_back(prototype, offset);
      }
    }
  }

  /**
   * @brief Zero out the ghost cell values on all the faces that have neighbors
   *
   * @tparam M the dimension of the faces
   * @param pinfo the PatchInfo for the patch
   * @param view the patch that is being zeroed
   */
  template<int M>
  void zeroGhostCellsOnAllFaces(const PatchInfo<D>& pinfo, PatchView<const double, D>& view) const
  {
    for (Face<D, M> f : Face<D, M>::getValues()) {
      if (pinfo.hasNbr(f)) {
        std::array<size_t, D - M> start;
        start.fill(0);
        std::array<size_t, D - M> end;
        end.fill(domain.getNumGhostCells() - 1);
        Loop::Nested<D - M>(start, end, [&](const std::array<size_t, D - M>& offset) {
          View<double, M + 1> this_ghost = view.getGhostSliceOn(f, offset);
          Loop::OverInteriorIndexes<M + 1>(
            this_ghost, [&](const std::array<int, M + 1>& coord) { this_ghost[coord] = 0; });
        });
      }
    }
  }

  /**
   * @brief Zero out the ghost cell values
   *
   * @param u the vector to zero out the ghost cells on
   */
  void zeroGhostCells(const Vector<D>& u) const
  {
    for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
      PatchView<const double, D> this_patch = u.getPatchView(pinfo.local_index);
      switch (fill_type) {
        case GhostFillingType::Corners:
          if constexpr (D >= 2) {
            zeroGhostCellsOnAllFaces<0>(pinfo, this_patch);
          }
        case GhostFillingType::Edges:
          if constexpr (D == 3) {
            zeroGhostCellsOnAllFaces<1>(pinfo, this_patch);
          }
        case GhostFillingType::Faces:
          zeroGhostCellsOnAllFaces<D - 1>(pinfo, this_patch);
          break;
        default:
          throw RuntimeError("Unsupported GhostFilling Type");
      }
    }
  }

public:
  /**
   * @brief Construct a new MPIGhostFiller object
   *
   * @param domain  the domain being used
   * @param fill_type  the number of side cases to address
   */
  MPIGhostFiller(const Domain<D>& domain, GhostFillingType fill_type)
    : domain(domain)
    , fill_type(fill_type)
  {
    std::map<int, RemoteCallSet> rank_to_remote_call_sets;

    switch (fill_type) {
      case GhostFillingType::Corners:
        if constexpr (D >= 2) {
          ghost_local_data_infos.template get<0>() = GhostViewInfo<0>(domain);
          enumerateCalls<0>(local_calls.template get<0>(), rank_to_remote_call_sets);
        }
      case GhostFillingType::Edges:
        if constexpr (D == 3) {
          ghost_local_data_infos.template get<1>() = GhostViewInfo<1>(domain);
          enumerateCalls<1>(local_calls.template get<1>(), rank_to_remote_call_sets);
        }
      case GhostFillingType::Faces:
        ghost_local_data_infos.template get<D - 1>() = GhostViewInfo<D - 1>(domain);
        enumerateCalls<D - 1>(local_calls.template get<D - 1>(), rank_to_remote_call_sets);
        break;
      default:
        throw RuntimeError("Unsupported GhostFilling Type");
    }

    remote_call_sets.reserve(rank_to_remote_call_sets.size());
    for (const auto& pair : rank_to_remote_call_sets) {
      remote_call_sets.push_back(pair.second);
    }
  }
  /**
   * @brief Fill the ghost cells for the neighboring patch
   *
   * @param pinfo the patch that ghost cells are being filled from
   * @param local_vew the view for patch that ghost cells are being filled from
   * @param nbr_view the view for the neighboring patch, where ghost cells are being
   * filled.
   * @param side the side that the neighboring patch is on
   * @param nbr_type the type of neighbor
   * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is
   * coarser
   */
  virtual void fillGhostCellsForNbrPatch(const PatchInfo<D>& pinfo,
                                         const PatchView<const double, D>& local_view,
                                         const PatchView<const double, D>& nbr_view,
                                         Side<D> side,
                                         NbrType nbr_type,
                                         Orthant<D - 1> orthant_on_coarse) const = 0;
  /**
   * @brief Fill the edge ghost cells for the neighboring patch
   *
   * @param pinfo the patch that ghost cells are being filled from
   * @param local_view the view for patch that ghost cells are being filled from
   * @param nbr_view  the view for the neighboring patch, where ghost cells are being
   * filled.
   * @param edge the edge that the neighboring patch is on
   * @param nbr_type the type of neighbor
   * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is
   * coarser
   */
  virtual void fillGhostCellsForEdgeNbrPatch(const PatchInfo<D>& pinfo,
                                             const PatchView<const double, D>& local_view,
                                             const PatchView<const double, D>& nbr_view,
                                             Edge edge,
                                             NbrType nbr_type,
                                             Orthant<1> orthant_on_coarse) const = 0;
  /**
   * @brief Fill the corner ghost cells for the neighboring patch
   *
   * @param pinfo the patch that ghost cells are being filled from
   * @param local_view the view for patch that ghost cells are being filled from
   * @param nbr_view  the view for the neighboring patch, where ghost cells are being
   * filled.
   * @param corner the edge that the neighboring patch is on
   * @param nbr_type the type of neighbor
   */
  virtual void fillGhostCellsForCornerNbrPatch(const PatchInfo<D>& pinfo,
                                               const PatchView<const double, D>& local_view,
                                               const PatchView<const double, D>& nbr_view,
                                               Corner<D> corner,
                                               NbrType nbr_type) const = 0;

  /**
   * @brief Perform any on this patches ghost cells.
   *
   * This may be necessary on some schemes because it needs data from the patch itself, not just
   * the neighboring patch
   *
   * @param pinfo the patch
   * @param view the view for the patch
   */
  virtual void fillGhostCellsForLocalPatch(const PatchInfo<D>& pinfo,
                                           const PatchView<const double, D>& view) const = 0;

  /**
   * @brief Fill ghost cells on a vector
   *
   * @param u  the vector
   */
  void fillGhost(const Vector<D>& u) const override
  {
    if constexpr (ENABLE_DEBUG) {
      if (u.getNumLocalPatches() != domain.getNumLocalPatches()) {
        throw RuntimeError("u vector is incorrect length. Expected Lenght of " +
                           std::to_string(domain.getNumLocalPatches()) + " but vector was length " +
                           std::to_string(u.getNumLocalPatches()));
      }
    }
    // zero out ghost cells
    zeroGhostCells(u);

    // allocate recv buffers and post recvs
    std::vector<std::vector<double>> recv_buffers(remote_call_sets.size());
    for (size_t i = 0; i < remote_call_sets.size(); i++) {
      recv_buffers[i].resize(remote_call_sets[i].recv_buffer_length * u.getNumComponents());
    }
    std::vector<MPI_Request> recv_requests = postRecvs(recv_buffers);

    // allocate send buffers
    std::vector<std::vector<double>> out_buffers(remote_call_sets.size());
    for (size_t i = 0; i < remote_call_sets.size(); i++) {
      out_buffers[i].resize(remote_call_sets[i].send_buffer_length * u.getNumComponents());
    }
    std::vector<MPI_Request> send_requests = postSends(out_buffers, u);

    // perform local operations
    for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
      PatchView<const double, D> view = u.getPatchView(pinfo.local_index);
      fillGhostCellsForLocalPatch(pinfo, view);
    }
    processLocalFills<D - 1>(u);

    processRecvs(recv_requests, recv_buffers, u);

    // wait for sends for finish
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUS_IGNORE);
  }

  /**
   * @brief Get the ghost filling type
   *
   * @return GhostFillingType the tyep
   */
  GhostFillingType getFillType() const { return fill_type; }

  /**
   * @brief Get the domain that is being filled for
   *
   * @return std::shared_ptr<const Domain<D>>  the domain
   */
  const Domain<D>& getDomain() const { return domain; }
};
extern template class MPIGhostFiller<1>;
extern template class MPIGhostFiller<2>;
extern template class MPIGhostFiller<3>;
} // namespace ThunderEgg
#endif
