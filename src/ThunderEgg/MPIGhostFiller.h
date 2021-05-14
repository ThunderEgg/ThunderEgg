/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
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

#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/GhostFillingType.h>
#include <mpi.h>
namespace ThunderEgg
{
/**
 * @brief Parallell ghostfiller implimented with MPI
 *
 * There are two private functions that have to be overridden derived classes.
 * fillGhostCellsForNbrPatch, and fillGhostCellsForLocalPatch
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class MPIGhostFiller : public GhostFiller<D>
{
	private:
	template <class Face> struct RemoteCallPrototype {
		/**
		 * @brief The id of the neighboring patch
		 */
		int id;
		/**
		 * @brief The face of the local patch that the neighboring patch is on
		 */
		Face face;
		/**
		 * @brief The type of neighbor
		 */
		NbrType nbr_type;
		/**
		 * @brief The orthant on the on the coarse patch, or null
		 */
		Orthant<Face::dimensionality> orthant;
		/**
		 * @brief The local index of the local patch the we are ghost filling from
		 */
		int local_index;
		RemoteCallPrototype(int id, Face face, NbrType nbr_type, Orthant<Face::dimensionality> orthant, int local_index)
		: id(id),
		  face(face),
		  nbr_type(nbr_type),
		  orthant(orthant),
		  local_index(local_index)
		{
		}
		bool operator<(const RemoteCallPrototype &other) const
		{
			return std::forward_as_tuple(id, face.opposite(), nbr_type, orthant, local_index)
			       < std::forward_as_tuple(other.id, other.face.opposite(), other.nbr_type, other.orthant, other.local_index);
		}
	};

	template <class Face> struct IncomingGhostPrototype {
		/**
		 * @brief The id of the patch that is being filled
		 */
		int id;
		/**
		 * @brief The side of the patch that ghosts are on
		 */
		Face face;
		/**
		 * @brief the local index of the patch being filled
		 */
		int local_index;
		IncomingGhostPrototype(int id, Face face, int local_index) : id(id), face(face), local_index(local_index) {}
		bool operator<(const IncomingGhostPrototype &other) const
		{
			return std::forward_as_tuple(id, face, local_index) < std::forward_as_tuple(other.id, other.face, other.local_index);
		}
	};

	/**
	 * @brief Struct for remote call
	 */
	template <class Face> struct RemoteCall {
		/**
		 * @brief The side of the patch that the neighboring patch is on
		 */
		Face face;
		/**
		 * @brief The type of neighbor
		 */
		NbrType nbr_type;
		/**
		 * @brief The orthant on the coarse patch, or null
		 */
		Orthant<Face::dimensionality> orthant;
		/**
		 * @brief The local index of the patch that we are filling from
		 *
		 */
		int local_index;
		/**
		 * @brief offset in buffer
		 */
		size_t offset;
		RemoteCall(const RemoteCallPrototype<Face> &prototype, size_t offset)
		: face(prototype.face),
		  nbr_type(prototype.nbr_type),
		  orthant(prototype.orthant),
		  local_index(prototype.local_index),
		  offset(offset)
		{
		}
	};

	template <class Face> struct IncomingGhost {
		/**
		 * @brief The side of the patch that the ghosts are on
		 */
		Face face;
		/**
		 * @brief The local index of the patch that is being filled
		 */
		int local_index;
		/**
		 * @brief The offset in the buffer
		 */
		size_t offset;
		IncomingGhost(const IncomingGhostPrototype<Face> &prototype, size_t offset)
		: face(prototype.face),
		  local_index(prototype.local_index),
		  offset(offset)
		{
		}
	};

	/**
	 * @brief Struct for local call
	 */
	template <class Face> struct LocalCall {
		/**
		 * @brief The side that the neighbor is on
		 */
		Face face;
		/**
		 * @brief The neighbor
		 *
		 */
		NbrType nbr_type;
		/**
		 * @brief The orthant on the coarse patch
		 */
		Orthant<Face::dimensionality> orthant;
		/**
		 * @brief The local index of the patch that is being filled from
		 */
		int local_index;
		/**
		 * @brief The local index of the patch that is being filled
		 */
		int nbr_local_index;
		LocalCall(Face face, NbrType nbr_type, Orthant<Face::dimensionality> orthant, int local_index, size_t nbr_local_index)
		: face(face),
		  nbr_type(nbr_type),
		  orthant(orthant),
		  local_index(local_index),
		  nbr_local_index(nbr_local_index)
		{
		}
	};

	struct RemoteCallSet {
		/**
		 * @brief The rank that we are communicating with
		 */
		int rank;
		/**
		 * @brief the length of the send buffer;
		 */
		size_t send_buffer_length = 0;
		/**
		 * @brief A vector of RemoteCall deques, one sperate deque for each rank
		 */
		std::deque<RemoteCall<Side<D>>> remote_calls;
		std::deque<RemoteCall<Edge>>    edge_remote_calls;
		/**
		 * @brief the length of the recv buffer;
		 */
		size_t recv_buffer_length = 0;
		/**
		 * @brief deque for incoming ghost calls
		 */
		std::deque<IncomingGhost<Side<D>>> incoming_ghosts;
		std::deque<IncomingGhost<Edge>>    edge_incoming_ghosts;

		explicit RemoteCallSet(int rank) : rank(rank) {}
	};

	std::vector<RemoteCallSet> remote_call_sets;
	/**
	 * @brief deque of local calls to be made
	 */
	std::deque<LocalCall<Side<D>>> local_calls;
	std::deque<LocalCall<Edge>>    edge_local_calls;

	/**
	 * @brief Get the LocalData object for the buffer
	 *
	 * @param buffer_ptr pointer to the ghost cells position in the buffer
	 * @param side  the side that the ghost cells are on
	 * @param component_index  the component index
	 * @return LocalData<D> the LocalData object
	 */
	LocalData<D> getLocalDataForBuffer(double *buffer_ptr, Side<D> side, int component_index) const
	{
		auto ns              = domain->getNs();
		int  num_ghost_cells = domain->getNumGhostCells();
		// determine striding
		std::array<int, D> strides;
		strides[0] = 1;
		for (size_t i = 1; i < D; i++) {
			if (i - 1 == side.getAxisIndex()) {
				strides[i] = num_ghost_cells * strides[i - 1];
			} else {
				strides[i] = ns[i - 1] * strides[i - 1];
			}
		}
		int size = D - 1 == side.getAxisIndex() ? (num_ghost_cells * strides[D - 1]) : (ns[D - 1] * strides[D - 1]);
		// transform buffer ptr so that it points to first non-ghost cell
		double *transformed_buffer_ptr = buffer_ptr + size * component_index;
		if (side.isLowerOnAxis()) {
			transformed_buffer_ptr -= (-num_ghost_cells) * strides[side.getAxisIndex()];
		} else {
			transformed_buffer_ptr -= ns[side.getAxisIndex()] * strides[side.getAxisIndex()];
		}

		LocalData<D> buffer_data(transformed_buffer_ptr, strides, ns, num_ghost_cells);
		return buffer_data;
	}
	/**
	 * @brief Post send requests
	 *
	 * @param buffers the allocated buffers
	 * @return std::vector<MPI_Request> the send requests
	 */
	std::vector<MPI_Request> postRecvs(std::vector<std::vector<double>> &buffers) const
	{
		std::vector<MPI_Request> recv_requests(buffers.size());
		for (size_t i = 0; i < recv_requests.size(); i++) {
			MPI_Irecv(buffers[i].data(), buffers[i].size(), MPI_DOUBLE, remote_call_sets[i].rank, 0, MPI_COMM_WORLD, &recv_requests[i]);
		}
		return recv_requests;
	}
	/**
	 * @brief process recv requests as they are ready
	 *
	 * @param requests the recv requests
	 * @param buffers the recv buffers
	 * @param u the vector to fill ghost values in
	 */
	void processRecvs(std::vector<MPI_Request> &requests, std::vector<std::vector<double>> &buffers, std::shared_ptr<const Vector<D>> u) const
	{
		size_t num_requests = requests.size();
		for (size_t i = 0; i < num_requests; i++) {
			int finished_index;
			MPI_Waitany(requests.size(), requests.data(), &finished_index, MPI_STATUS_IGNORE);
			for (const IncomingGhost<Side<D>> &incoming_ghost : remote_call_sets[finished_index].incoming_ghosts) {
				int     local_index   = incoming_ghost.local_index;
				Side<D> side          = incoming_ghost.face;
				size_t  buffer_offset = incoming_ghost.offset;

				for (int c = 0; c < u->getNumComponents(); c++) {
					const LocalData<D> local_data  = u->getLocalData(c, local_index);
					double *           buffer_ptr  = buffers[finished_index].data() + buffer_offset * u->getNumComponents();
					LocalData<D>       buffer_data = getLocalDataForBuffer(buffer_ptr, side, c);
					for (int ig = 0; ig < domain->getNumGhostCells(); ig++) {
						LocalData<D - 1> local_slice  = local_data.getGhostSliceOnSide(side, ig + 1);
						LocalData<D - 1> buffer_slice = buffer_data.getGhostSliceOnSide(side, ig + 1);
						nested_loop<D - 1>(local_slice.getStart(), local_slice.getEnd(), [&](const std::array<int, D - 1> &coord) {
							local_slice[coord] += buffer_slice[coord];
						});
					}
				}
			}
		}
	}
	/**
	 * @brief fill buffers and post send requests
	 *
	 * @param buffers the allocated buffers
	 * @param u the vector to fill buffers from
	 * @return std::vector<MPI_Request> the requests
	 */
	std::vector<MPI_Request> postSends(std::vector<std::vector<double>> &buffers, std::shared_ptr<const Vector<D>> u) const
	{
		std::vector<MPI_Request> send_requests(remote_call_sets.size());
		for (size_t i = 0; i < remote_call_sets.size(); i++) {
			for (const RemoteCall<Side<D>> &call : remote_call_sets[i].remote_calls) {
				const PatchInfo<D> &pinfo       = domain->getPatchInfoVector()[call.local_index];
				Side<D>             side        = call.face.opposite();
				auto                local_datas = u->getLocalDatas(call.local_index);
				double *            buffer_ptr  = buffers[i].data() + call.offset * u->getNumComponents();

				// create LocalData objects for the buffer
				std::vector<LocalData<D>> buffer_datas;
				buffer_datas.reserve(u->getNumComponents());
				for (int c = 0; c < u->getNumComponents(); c++) {
					buffer_datas.push_back(getLocalDataForBuffer(buffer_ptr, side, c));
				}

				// make the call
				fillGhostCellsForNbrPatch(pinfo, local_datas, buffer_datas, call.face, call.nbr_type, call.orthant);
			}
			MPI_Isend(buffers[i].data(), buffers[i].size(), MPI_DOUBLE, remote_call_sets[i].rank, 0, MPI_COMM_WORLD, &send_requests[i]);
		}
		return send_requests;
	}

	template <class Face>
	static void addNormalNbrCalls(const PatchInfo<D> &                                   pinfo,
	                              Face                                                   f,
	                              std::deque<LocalCall<Face>> &                          my_local_calls,
	                              std::map<int, std::set<RemoteCallPrototype<Face>>> &   rank_to_remote_call_prototypes,
	                              std::map<int, std::set<IncomingGhostPrototype<Face>>> &rank_to_incoming_ghost_prototypes)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto nbrinfo = pinfo.getNormalNbrInfo(f);
		if (nbrinfo.rank == rank) {
			my_local_calls.emplace_back(f, NbrType::Normal, Orthant<Face::dimensionality>::null(), pinfo.local_index, nbrinfo.local_index);
		} else {
			rank_to_remote_call_prototypes[nbrinfo.rank].emplace(
			nbrinfo.id, f, NbrType::Normal, Orthant<Face::dimensionality>::null(), pinfo.local_index);
			rank_to_incoming_ghost_prototypes[nbrinfo.rank].emplace(pinfo.id, f, pinfo.local_index);
		}
	}

	template <class Face>
	static void addFineNbrCalls(const PatchInfo<D> &                                   pinfo,
	                            Face                                                   f,
	                            std::deque<LocalCall<Face>> &                          my_local_calls,
	                            std::map<int, std::set<RemoteCallPrototype<Face>>> &   rank_to_remote_call_prototypes,
	                            std::map<int, std::set<IncomingGhostPrototype<Face>>> &rank_to_incoming_ghost_prototypes)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto nbrinfo = pinfo.getFineNbrInfo(f);
		for (size_t i = 0; i < Orthant<Face::dimensionality>::num_orthants; i++) {
			if (nbrinfo.ranks[i] == rank) {
				my_local_calls.emplace_back(f, NbrType::Fine, Orthant<Face::dimensionality>(i), pinfo.local_index, nbrinfo.local_indexes[i]);
			} else {
				rank_to_remote_call_prototypes[nbrinfo.ranks[i]].emplace(
				nbrinfo.ids[i], f, NbrType::Fine, Orthant<Face::dimensionality>(i), pinfo.local_index);
				rank_to_incoming_ghost_prototypes[nbrinfo.ranks[i]].emplace(pinfo.id, f, pinfo.local_index);
			}
		}
	}

	template <class Face>
	static void addCoarseNbrCalls(const PatchInfo<D> &                                   pinfo,
	                              Face                                                   f,
	                              std::deque<LocalCall<Face>> &                          my_local_calls,
	                              std::map<int, std::set<RemoteCallPrototype<Face>>> &   rank_to_remote_call_prototypes,
	                              std::map<int, std::set<IncomingGhostPrototype<Face>>> &rank_to_incoming_ghost_prototypes)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto nbrinfo = pinfo.getCoarseNbrInfo(f);
		if (nbrinfo.rank == rank) {
			my_local_calls.emplace_back(f, NbrType::Coarse, nbrinfo.orth_on_coarse, pinfo.local_index, nbrinfo.local_index);
		} else {
			rank_to_remote_call_prototypes[nbrinfo.rank].emplace(nbrinfo.id, f, NbrType::Coarse, nbrinfo.orth_on_coarse, pinfo.local_index);
			rank_to_incoming_ghost_prototypes[nbrinfo.rank].emplace(pinfo.id, f, pinfo.local_index);
		}
	}
	template <class Face> std::vector<size_t> getGhostLengths() const
	{
		std::vector<size_t> ghost_lengths(12);
		if constexpr (std::is_same<Face, Side<D>>::value) {
			for (size_t axis = 0; axis < D; axis++) {
				size_t length = 1;
				for (size_t i = 0; i < D; i++) {
					if (i == axis) {
						length *= domain->getNumGhostCells();
					} else {
						length *= domain->getNs()[i];
					}
				}
				ghost_lengths[axis] = length;
			}
		} else if constexpr (std::is_same<Face, Edge>::value) {
			for (size_t axis = 0; axis < D; axis++) {
				size_t length = 1;
				for (size_t i = 0; i < D; i++) {
					if (i == axis) {
						length *= domain->getNs()[i];
					} else {
						length *= domain->getNumGhostCells();
					}
				}
				ghost_lengths[axis] = length;
			}
		} else {
			static_assert(std::is_same_v<Face *, void>, "Unsupported type");
		}
		return ghost_lengths;
	}
	template <class Face> size_t getGhostLength(const std::vector<size_t> &ghost_lengths, Face face) const
	{
		size_t length = 0;
		if constexpr (std::is_same<Face, Side<D>>::value) {
			if constexpr (std::is_same<Face, Side<D>>::value || std::is_same<Face, Edge>::value) {
				length = ghost_lengths[face.getAxisIndex()];
			} else {
				static_assert(std::is_same_v<Face *, void>, "Unsupported type");
			}
		}
		return length;
	}
	template <class Face>
	void enumerateCalls(std::deque<LocalCall<Face>> &my_local_calls, std::map<int, RemoteCallSet> &rank_to_remote_call_sets) const
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::map<int, std::set<RemoteCallPrototype<Face>>>    rank_to_remote_call_prototypes;
		std::map<int, std::set<IncomingGhostPrototype<Face>>> rank_to_incoming_ghost_prototypes;

		for (const PatchInfo<D> &pinfo : domain->getPatchInfoVector()) {
			for (Face f : Face::getValues()) {
				if (pinfo.hasNbr(f)) {
					switch (pinfo.getNbrType(f)) {
						case NbrType::Normal:
							addNormalNbrCalls(pinfo, f, my_local_calls, rank_to_remote_call_prototypes, rank_to_incoming_ghost_prototypes);
							break;
						case NbrType::Fine:
							addFineNbrCalls(pinfo, f, my_local_calls, rank_to_remote_call_prototypes, rank_to_incoming_ghost_prototypes);
							break;
						case NbrType::Coarse:
							addCoarseNbrCalls(pinfo, f, my_local_calls, rank_to_remote_call_prototypes, rank_to_incoming_ghost_prototypes);
							break;
						default:
							throw RuntimeError("Unsupported NbrType");
					}
				}
			}
		}
		std::vector<size_t> ghost_lengths = getGhostLengths<Face>();
		for (const auto &pair : rank_to_remote_call_prototypes) {
			int rank = pair.first;
			rank_to_remote_call_sets.emplace(rank, rank);
			RemoteCallSet &       remote_call_set = rank_to_remote_call_sets.at(rank);
			std::tuple<int, Face> prev_id_side;
			for (const RemoteCallPrototype<Face> &call : pair.second) {
				size_t offset = remote_call_set.send_buffer_length;

				// calculate length in buffer need for ghost cells
				size_t length = getGhostLength(ghost_lengths, call.face);
				// add length to buffer length
				// if its the same side of the patch, resuse the previous buffer space
				if (std::make_tuple(call.id, call.face) == prev_id_side) {
					offset -= length;
				} else {
					remote_call_set.send_buffer_length += length;
				}
				if constexpr (std::is_same<Face, Side<D>>::value) {
					remote_call_set.remote_calls.emplace_back(call, offset);
				} else if constexpr (std::is_same<Face, Edge>::value) {
					remote_call_set.edge_remote_calls.emplace_back(call, offset);
				}
				prev_id_side = std::make_tuple(call.id, call.face);
			}
		}
		for (const auto &pair : rank_to_incoming_ghost_prototypes) {
			RemoteCallSet &remote_call_set = rank_to_remote_call_sets.at(pair.first);
			for (const IncomingGhostPrototype<Face> &prototype : pair.second) {
				// add length for ghosts to buffer length
				size_t length = getGhostLength(ghost_lengths, prototype.face);
				size_t offset = remote_call_set.recv_buffer_length;
				remote_call_set.recv_buffer_length += length;

				// add ghost to incoming ghosts
				if constexpr (std::is_same<Face, Side<D>>::value) {
					remote_call_set.incoming_ghosts.emplace_back(prototype, offset);
				} else if constexpr (std::is_same<Face, Edge>::value) {
					remote_call_set.edge_incoming_ghosts.emplace_back(prototype, offset);
				}
			}
		}
	}

	protected:
	/**
	 * @brief The domain that this ghostfiller operates on
	 */
	std::shared_ptr<const Domain<D>> domain;
	/**
	 * @brief the fill type to use
	 */
	GhostFillingType fill_type;

	public:
	/**
	 * @brief Construct a new MPIGhostFiller object
	 *
	 * @param domain  the domain being used
	 * @param fill_type  the number of side cases to address
	 */
	MPIGhostFiller(std::shared_ptr<const Domain<D>> domain, GhostFillingType fill_type) : domain(domain), fill_type(fill_type)
	{
		std::map<int, RemoteCallSet> rank_to_remote_call_sets;

		switch (fill_type) {
			case GhostFillingType::Corners:
			case GhostFillingType::Edges:
				if constexpr (D == 3) {
					enumerateCalls<Edge>(edge_local_calls, rank_to_remote_call_sets);
				}
			case GhostFillingType::Faces:
				enumerateCalls<Side<D>>(local_calls, rank_to_remote_call_sets);
				break;
			default:
				throw RuntimeError("Unsupported GhostFilling Type");
		}

		remote_call_sets.reserve(rank_to_remote_call_sets.size());
		for (const auto &pair : rank_to_remote_call_sets) {
			remote_call_sets.push_back(pair.second);
		}
	}
	/**
	 * @brief Fill the ghost cells for the neighboring patch
	 *
	 * @param pinfo the patch that ghost cells are being filled from
	 * @param local_datas the local data for patch that ghost cells are being filled from
	 * @param nbr_datas  the local data for the neighboring patch, where ghost cells are being
	 * filled.
	 * @param side the side that the neighboring patch is on
	 * @param nbr_type the type of neighbor
	 * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is coarser
	 */
	virtual void fillGhostCellsForNbrPatch(const PatchInfo<D> &             pinfo,
	                                       const std::vector<LocalData<D>> &local_datas,
	                                       std::vector<LocalData<D>> &      nbr_datas,
	                                       Side<D>                          side,
	                                       NbrType                          nbr_type,
	                                       Orthant<D - 1>                   orthant_on_coarse) const = 0;
	/**
	 * @brief Fill the edge ghost cells for the neighboring patch
	 *
	 * @param pinfo the patch that ghost cells are being filled from
	 * @param local_datas the local data for patch that ghost cells are being filled from
	 * @param nbr_datas  the local data for the neighboring patch, where ghost cells are being
	 * filled.
	 * @param edge the edge that the neighboring patch is on
	 * @param nbr_type the type of neighbor
	 * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is coarser
	 */
	virtual void fillGhostCellsForEdgeNbrPatch(const PatchInfo<D> &             pinfo,
	                                           const std::vector<LocalData<D>> &local_datas,
	                                           std::vector<LocalData<D>> &      nbr_datas,
	                                           Edge                             edge,
	                                           NbrType                          nbr_type,
	                                           Orthant<1>                       orthant_on_coarse) const = 0;
	/**
	 * @brief Fill the corner ghost cells for the neighboring patch
	 *
	 * @param pinfo the patch that ghost cells are being filled from
	 * @param local_datas the local data for patch that ghost cells are being filled from
	 * @param nbr_datas  the local data for the neighboring patch, where ghost cells are being
	 * filled.
	 * @param corner the edge that the neighboring patch is on
	 * @param nbr_type the type of neighbor
	 */
	virtual void fillGhostCellsForCornerNbrPatch(const PatchInfo<D> &             pinfo,
	                                             const std::vector<LocalData<D>> &local_datas,
	                                             std::vector<LocalData<D>> &      nbr_datas,
	                                             Corner<D>                        corner,
	                                             NbrType                          nbr_type) const = 0;

	/**
	 * @brief Perform any on this patches ghost cells.
	 *
	 * This may be necessary on some schemes because it needs data from the patch itself, not just
	 * the neighboring patch
	 *
	 * @param pinfo the patch
	 * @param local_datas the LocalData for the patch
	 */
	virtual void fillGhostCellsForLocalPatch(const PatchInfo<D> &pinfo, std::vector<LocalData<D>> &local_datas) const = 0;

	/**
	 * @brief Fill ghost cells on a vector
	 *
	 * @param u  the vector
	 */
	void fillGhost(std::shared_ptr<const Vector<D>> u) const
	{
		// zero out ghost cells
		for (const PatchInfo<D> &pinfo : domain->getPatchInfoVector()) {
			for (auto &this_patch : u->getLocalDatas(pinfo.local_index)) {
				for (Side<D> s : Side<D>::getValues()) {
					if (pinfo.hasNbr(s)) {
						for (int i = 0; i < pinfo.num_ghost_cells; i++) {
							auto this_ghost = this_patch.getGhostSliceOnSide(s, i + 1);
							nested_loop<D - 1>(
							this_ghost.getStart(), this_ghost.getEnd(), [&](const std::array<int, D - 1> &coord) { this_ghost[coord] = 0; });
						}
					}
				}
				if constexpr (D == 3) {
					for (Edge e : Edge::getValues()) {
						if (pinfo.hasNbr(e)) {
							for (int j = 0; j < pinfo.num_ghost_cells; j++) {
								for (int i = 0; i < pinfo.num_ghost_cells; i++) {
									auto this_ghost = this_patch.getSliceOnEdge(e, {-1 - i, -1 - j});
									nested_loop<1>(
									this_ghost.getStart(), this_ghost.getEnd(), [&](const std::array<int, 1> &coord) { this_ghost[coord] = 0; });
								}
							}
						}
					}
				}
			}
		}

		// allocate recv buffers and post recvs
		std::vector<std::vector<double>> recv_buffers(remote_call_sets.size());
		for (size_t i = 0; i < remote_call_sets.size(); i++) {
			recv_buffers[i].resize(remote_call_sets[i].recv_buffer_length * u->getNumComponents());
		}
		std::vector<MPI_Request> recv_requests = postRecvs(recv_buffers);

		// allocate send buffers
		std::vector<std::vector<double>> out_buffers(remote_call_sets.size());
		for (size_t i = 0; i < remote_call_sets.size(); i++) {
			out_buffers[i].resize(remote_call_sets[i].send_buffer_length * u->getNumComponents());
		}
		std::vector<MPI_Request> send_requests = postSends(out_buffers, u);

		// perform local operations
		for (const PatchInfo<D> &pinfo : domain->getPatchInfoVector()) {
			auto datas = u->getLocalDatas(pinfo.local_index);
			fillGhostCellsForLocalPatch(pinfo, datas);
		}
		for (const LocalCall<Side<D>> &call : local_calls) {
			const PatchInfo<D> &pinfo       = domain->getPatchInfoVector()[call.local_index];
			auto                local_datas = u->getLocalDatas(call.local_index);
			auto                nbr_datas   = u->getLocalDatas(call.nbr_local_index);
			fillGhostCellsForNbrPatch(pinfo, local_datas, nbr_datas, call.face, call.nbr_type, call.orthant);
		}
		if constexpr (D == 3) {
			for (const LocalCall<Edge> &call : edge_local_calls) {
				const PatchInfo<D> &pinfo       = domain->getPatchInfoVector()[call.local_index];
				auto                local_datas = u->getLocalDatas(call.local_index);
				auto                nbr_datas   = u->getLocalDatas(call.nbr_local_index);
				fillGhostCellsForEdgeNbrPatch(pinfo, local_datas, nbr_datas, call.face, call.nbr_type, call.orthant);
			}
		}

		processRecvs(recv_requests, recv_buffers, u);

		// wait for sends for finish
		MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUS_IGNORE);
	}
};
extern template class MPIGhostFiller<1>;
extern template class MPIGhostFiller<2>;
extern template class MPIGhostFiller<3>;
} // namespace ThunderEgg
#endif
