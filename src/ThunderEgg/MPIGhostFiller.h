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
	template <int N, template <int> class T> class DimensionalArray
	{
		private:
		T<N>                       t;
		DimensionalArray<N - 1, T> next;

		public:
		template <int I> T<I> &get()
		{
			static_assert(I < N, "invalid index value");
			return next.template get<I>();
		}
		template <> T<N> &get()
		{
			return t;
		}
		template <int I> const T<I> &get() const
		{
			static_assert(I < N, "invalid index value");
			return next.template get<I>();
		}
		template <> const T<N> &get() const
		{
			return t;
		}
	};
	template <template <int> class T> class DimensionalArray<0, T>
	{
		private:
		T<0> t;

		public:
		template <int I> T<I> &get()
		{
			static_assert(I == 0, "invalid index value");
		}
		template <> T<0> &get()
		{
			return t;
		}
		template <int I> const T<I> &get() const
		{
			static_assert(I == 0, "invalid index value");
		}
		template <> const T<0> &get() const
		{
			return t;
		}
	};
	template <int M> struct RemoteCallPrototype {
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
		RemoteCallPrototype(int id, Face<D, M> face, NbrType nbr_type, Orthant<M> orthant, int local_index)
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

	template <int M> struct IncomingGhostPrototype {
		/**
		 * @brief The id of the patch that is being filled
		 */
		int id;
		/**
		 * @brief The side of the patch that ghosts are on
		 */
		Face<D, M> face;
		/**
		 * @brief the local index of the patch being filled
		 */
		int local_index;
		IncomingGhostPrototype(int id, Face<D, M> face, int local_index) : id(id), face(face), local_index(local_index) {}
		bool operator<(const IncomingGhostPrototype &other) const
		{
			return std::forward_as_tuple(id, face, local_index) < std::forward_as_tuple(other.id, other.face, other.local_index);
		}
	};

	/**
	 * @brief Struct for remote call
	 */
	template <int M> struct RemoteCall {
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
		 * @brief offset in buffer
		 */
		size_t offset;
		RemoteCall(const RemoteCallPrototype<M> &prototype, size_t offset)
		: face(prototype.face),
		  nbr_type(prototype.nbr_type),
		  orthant(prototype.orthant),
		  local_index(prototype.local_index),
		  offset(offset)
		{
		}
	};

	template <int M> struct IncomingGhost {
		/**
		 * @brief The side of the patch that the ghosts are on
		 */
		Face<D, M> face;
		/**
		 * @brief The local index of the patch that is being filled
		 */
		int local_index;
		/**
		 * @brief The offset in the buffer
		 */
		size_t offset;
		IncomingGhost(const IncomingGhostPrototype<M> &prototype, size_t offset)
		: face(prototype.face),
		  local_index(prototype.local_index),
		  offset(offset)
		{
		}
	};

	/**
	 * @brief Struct for local call
	 */
	template <int M> struct LocalCall {
		/**
		 * @brief The side that the neighbor is on
		 */
		Face<D, M> face;
		/**
		 * @brief The neighbor
		 *
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
		LocalCall(Face<D, M> face, NbrType nbr_type, Orthant<M> orthant, int local_index, size_t nbr_local_index)
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
		size_t send_buffer_length              = 0;
		template <int M> using RemoteCallDeque = std::deque<RemoteCall<M>>;
		/**
		 * @brief A vector of RemoteCall deques, one sperate deque for each rank
		 */
		DimensionalArray<D - 1, RemoteCallDeque> remote_calls;
		/**
		 * @brief the length of the recv buffer;
		 */
		size_t recv_buffer_length                 = 0;
		template <int M> using IncomingGhostDeque = std::deque<IncomingGhost<M>>;
		/**
		 * @brief deque for incoming ghost calls
		 */
		DimensionalArray<D - 1, IncomingGhostDeque> incoming_ghosts;

		explicit RemoteCallSet(int rank) : rank(rank) {}
	};

	std::vector<RemoteCallSet> remote_call_sets;
	template <int M> using LocalCallDeque = std::deque<LocalCall<M>>;
	/**
	 * @brief deque of local calls to be made
	 */
	DimensionalArray<D - 1, LocalCallDeque> local_calls;

	template <int M> class GhostLocalDataInfo
	{
		private:
		std::array<std::array<int, D>, Face<D, M>::number_of> strides;
		std::array<int, Face<D, M>::number_of>                start_offsets;
		std::array<size_t, Face<D, M>::number_of>             sizes;

		public:
		GhostLocalDataInfo() = default;
		GhostLocalDataInfo(const Domain<D> &domain)
		{
			std::array<int, D> ns              = domain.getNs();
			int                num_ghost_cells = domain.getNumGhostCells();
			for (Face<D, M> face : Face<D, M>::getValues()) {
				std::array<int, D> ghost_ns = ns;
				for (Side<D> side : face.getSides()) {
					ghost_ns[side.getAxisIndex()] = num_ghost_cells;
				}
				strides[face.getIndex()][0] = 1;
				for (size_t i = 1; i < D; i++) {
					strides[face.getIndex()][i] = ghost_ns[i - 1] * strides[face.getIndex()][i - 1];
				}
				sizes[face.getIndex()] = ghost_ns[D - 1] * strides[face.getIndex()][D - 1];
				int offset             = 0;
				for (Side<D> side : face.getSides()) {
					if (side.isLowerOnAxis()) {
						offset += num_ghost_cells * strides[face.getIndex()][side.getAxisIndex()];
					} else {
						offset += -ns[side.getAxisIndex()] * strides[face.getIndex()][side.getAxisIndex()];
					}
				}
				start_offsets[face.getIndex()] = offset;
				/*
				for (int axis = 0; axis < D; axis++) {
				    if (sides[sides_index].getAxisIndex() == axis) {
				        if (sides[sides_index].isLowerOnAxis()) {
				            new_data += offset[sides_index] * strides[axis];
				        } else {
				            new_data += (lengths[axis] - 1 - offset[sides_index]) * strides[axis];
				        }
				        sides_index++;
				    } else {
				        new_lengths[lengths_index] = lengths[axis];
				        new_strides[lengths_index] = strides[axis];
				        lengths_index++;
				    }
				}
				*/
			}
		}
		size_t getSize(Face<D, M> face) const
		{
			return sizes[face.getIndex()];
		}
		const std::array<int, D> &getStrides(Face<D, M> face) const
		{
			return strides[face.getIndex()];
		}
		size_t getStartOffset(Face<D, M> face) const
		{
			return start_offsets[face.getIndex()];
		}
	};
	DimensionalArray<D - 1, GhostLocalDataInfo> ghost_local_data_infos;

	/**
	 * @brief Get the LocalData object for the buffer
	 *
	 * @param buffer_ptr pointer to the ghost cells position in the buffer
	 * @param side  the side that the ghost cells are on
	 * @param component_index  the component index
	 * @return LocalData<D> the LocalData object
	 */
	template <int M> LocalData<D> getLocalDataForBuffer(double *buffer_ptr, Face<D, M> face, int component_index) const
	{
		const GhostLocalDataInfo<M> &gld_info = ghost_local_data_infos.template get<M>();
		return LocalData<D>(buffer_ptr + gld_info.getSize(face) * component_index + gld_info.getStartOffset(face),
		                    gld_info.getStrides(face),
		                    domain->getNs(),
		                    domain->getNumGhostCells());
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
	template <int M> void addRecvBufferToGhost(const RemoteCallSet &remote_call_set, std::vector<double> &buffer, const Vector<D> &u) const
	{
		for (const IncomingGhost<M> &incoming_ghost : remote_call_set.incoming_ghosts.template get<M>()) {
			int        local_index   = incoming_ghost.local_index;
			Face<D, M> face          = incoming_ghost.face;
			size_t     buffer_offset = incoming_ghost.offset;

			for (int c = 0; c < u.getNumComponents(); c++) {
				const LocalData<D>        local_data  = u.getLocalData(c, local_index);
				double *                  buffer_ptr  = buffer.data() + buffer_offset * u.getNumComponents();
				LocalData<D>              buffer_data = getLocalDataForBuffer(buffer_ptr, face, c);
				std::array<size_t, D - M> start;
				start.fill(0);
				std::array<size_t, D - M> end;
				end.fill(domain->getNumGhostCells() - 1);
				nested_loop<D - M>(start, end, [&](const std::array<size_t, D - M> &offset) {
					LocalData<M> local_slice  = local_data.getGhostSliceOn(face, offset);
					LocalData<M> buffer_slice = buffer_data.getGhostSliceOn(face, offset);
					nested_loop<M>(local_slice.getStart(), local_slice.getEnd(), [&](const std::array<int, M> &coord) {
						local_slice[coord] += buffer_slice[coord];
					});
				});
			}
		}
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
			switch (fill_type) {
				case GhostFillingType::Corners:
					if constexpr (D >= 2) {
						addRecvBufferToGhost<0>(remote_call_sets[finished_index], buffers[finished_index], *u);
					}
				case GhostFillingType::Edges:
					if constexpr (D == 3) {
						addRecvBufferToGhost<1>(remote_call_sets[finished_index], buffers[finished_index], *u);
					}
				case GhostFillingType::Faces:
					addRecvBufferToGhost<D - 1>(remote_call_sets[finished_index], buffers[finished_index], *u);
					break;
				default:
					throw RuntimeError("Unsupported GhostFilling Type");
			}
		}
	}

	template <int M> void fillSendBuffer(const RemoteCallSet &remote_call_set, std::vector<double> &buffer, const Vector<D> &u) const
	{
		for (const RemoteCall<M> &call : remote_call_set.remote_calls.template get<M>()) {
			const PatchInfo<D> &pinfo       = domain->getPatchInfoVector()[call.local_index];
			Face<D, M>          face        = call.face.opposite();
			auto                local_datas = u.getLocalDatas(call.local_index);
			double *            buffer_ptr  = buffer.data() + call.offset * u.getNumComponents();

			// create LocalData objects for the buffer
			std::vector<LocalData<D>> buffer_datas;
			buffer_datas.reserve(u.getNumComponents());
			for (int c = 0; c < u.getNumComponents(); c++) {
				buffer_datas.push_back(getLocalDataForBuffer(buffer_ptr, face, c));
			}

			// make the call
			fillGhostCellsForNbrPatchPriv(pinfo, local_datas, buffer_datas, call.face, call.nbr_type, call.orthant);
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
			switch (fill_type) {
				case GhostFillingType::Corners:
					if constexpr (D >= 2) {
						fillSendBuffer<0>(remote_call_sets[i], buffers[i], *u);
					}
				case GhostFillingType::Edges:
					if constexpr (D == 3) {
						fillSendBuffer<1>(remote_call_sets[i], buffers[i], *u);
					}
				case GhostFillingType::Faces:
					fillSendBuffer<D - 1>(remote_call_sets[i], buffers[i], *u);
					break;
				default:
					throw RuntimeError("Unsupported GhostFilling Type");
			}
			MPI_Isend(buffers[i].data(), buffers[i].size(), MPI_DOUBLE, remote_call_sets[i].rank, 0, MPI_COMM_WORLD, &send_requests[i]);
		}
		return send_requests;
	}

	/**
	 * @brief Fill the ghost cells for the neighboring patch
	 *
	 * @param pinfo the patch that ghost cells are being filled from
	 * @param local_datas the local data for patch that ghost cells are being filled from
	 * @param nbr_datas  the local data for the neighboring patch, where ghost cells are being
	 * filled.
	 * @param face the face that the neighboring patch is on
	 * @param nbr_type the type of neighbor
	 * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is coarser
	 */
	template <int M>
	void fillGhostCellsForNbrPatchPriv(const PatchInfo<D> &             pinfo,
	                                   const std::vector<LocalData<D>> &local_datas,
	                                   std::vector<LocalData<D>> &      nbr_datas,
	                                   Face<D, M>                       face,
	                                   NbrType                          nbr_type,
	                                   Orthant<M>                       orthant_on_coarse) const
	{
		if constexpr (M == D - 1) {
			fillGhostCellsForNbrPatch(pinfo, local_datas, nbr_datas, face, nbr_type, orthant_on_coarse);
		} else if constexpr (M > 0) {
			fillGhostCellsForEdgeNbrPatch(pinfo, local_datas, nbr_datas, face, nbr_type, orthant_on_coarse);
		} else {
			fillGhostCellsForCornerNbrPatch(pinfo, local_datas, nbr_datas, face, nbr_type);
		}
	}

	template <int M> void processLocalFills(std::shared_ptr<const Vector<D>> u) const
	{
		for (const LocalCall<M> &call : local_calls.template get<M>()) {
			const PatchInfo<D> &pinfo       = domain->getPatchInfoVector()[call.local_index];
			auto                local_datas = u->getLocalDatas(call.local_index);
			auto                nbr_datas   = u->getLocalDatas(call.nbr_local_index);
			fillGhostCellsForNbrPatchPriv(pinfo, local_datas, nbr_datas, call.face, call.nbr_type, call.orthant);
		}
		if constexpr (M > 0) {
			processLocalFills<M - 1>(u);
		}
	}

	template <int M>
	static void addNormalNbrCalls(const PatchInfo<D> &                                pinfo,
	                              Face<D, M>                                          f,
	                              std::deque<LocalCall<M>> &                          my_local_calls,
	                              std::map<int, std::set<RemoteCallPrototype<M>>> &   rank_to_remote_call_prototypes,
	                              std::map<int, std::set<IncomingGhostPrototype<M>>> &rank_to_incoming_ghost_prototypes)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto nbrinfo = pinfo.getNormalNbrInfo(f);
		if (nbrinfo.rank == rank) {
			my_local_calls.emplace_back(f, NbrType::Normal, Orthant<M>::null(), pinfo.local_index, nbrinfo.local_index);
		} else {
			rank_to_remote_call_prototypes[nbrinfo.rank].emplace(nbrinfo.id, f, NbrType::Normal, Orthant<M>::null(), pinfo.local_index);
			rank_to_incoming_ghost_prototypes[nbrinfo.rank].emplace(pinfo.id, f, pinfo.local_index);
		}
	}

	template <int M>
	static void addFineNbrCalls(const PatchInfo<D> &                                pinfo,
	                            Face<D, M>                                          f,
	                            std::deque<LocalCall<M>> &                          my_local_calls,
	                            std::map<int, std::set<RemoteCallPrototype<M>>> &   rank_to_remote_call_prototypes,
	                            std::map<int, std::set<IncomingGhostPrototype<M>>> &rank_to_incoming_ghost_prototypes)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto nbrinfo = pinfo.getFineNbrInfo(f);
		for (size_t i = 0; i < Orthant<M>::num_orthants; i++) {
			if (nbrinfo.ranks[i] == rank) {
				my_local_calls.emplace_back(f, NbrType::Fine, Orthant<M>(i), pinfo.local_index, nbrinfo.local_indexes[i]);
			} else {
				rank_to_remote_call_prototypes[nbrinfo.ranks[i]].emplace(nbrinfo.ids[i], f, NbrType::Fine, Orthant<M>(i), pinfo.local_index);
				rank_to_incoming_ghost_prototypes[nbrinfo.ranks[i]].emplace(pinfo.id, f, pinfo.local_index);
			}
		}
	}

	template <int M>
	static void addCoarseNbrCalls(const PatchInfo<D> &                                pinfo,
	                              Face<D, M>                                          f,
	                              std::deque<LocalCall<M>> &                          my_local_calls,
	                              std::map<int, std::set<RemoteCallPrototype<M>>> &   rank_to_remote_call_prototypes,
	                              std::map<int, std::set<IncomingGhostPrototype<M>>> &rank_to_incoming_ghost_prototypes)
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

	template <int M> void enumerateCalls(std::deque<LocalCall<M>> &my_local_calls, std::map<int, RemoteCallSet> &rank_to_remote_call_sets) const
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::map<int, std::set<RemoteCallPrototype<M>>>    rank_to_remote_call_prototypes;
		std::map<int, std::set<IncomingGhostPrototype<M>>> rank_to_incoming_ghost_prototypes;

		for (const PatchInfo<D> &pinfo : domain->getPatchInfoVector()) {
			for (Face<D, M> f : Face<D, M>::getValues()) {
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
		const GhostLocalDataInfo<M> &ghost_local_data_info = ghost_local_data_infos.template get<M>();
		for (const auto &pair : rank_to_remote_call_prototypes) {
			int rank = pair.first;
			rank_to_remote_call_sets.emplace(rank, rank);
			RemoteCallSet &             remote_call_set = rank_to_remote_call_sets.at(rank);
			std::tuple<int, Face<D, M>> prev_id_side;
			for (const RemoteCallPrototype<M> &call : pair.second) {
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
		for (const auto &pair : rank_to_incoming_ghost_prototypes) {
			RemoteCallSet &remote_call_set = rank_to_remote_call_sets.at(pair.first);
			for (const IncomingGhostPrototype<M> &prototype : pair.second) {
				// add length for ghosts to buffer length
				size_t length = ghost_local_data_info.getSize(prototype.face);
				size_t offset = remote_call_set.recv_buffer_length;
				remote_call_set.recv_buffer_length += length;

				// add ghost to incoming ghosts
				remote_call_set.incoming_ghosts.template get<M>().emplace_back(prototype, offset);
			}
		}
	}

	template <int M> void zeroGhostCellsOnAllFaces(const PatchInfo<D> &pinfo, const LocalData<D> &data) const
	{
		for (Face<D, M> f : Face<D, M>::getValues()) {
			if (pinfo.hasNbr(f)) {
				std::array<size_t, D - M> start;
				start.fill(0);
				std::array<size_t, D - M> end;
				end.fill(domain->getNumGhostCells() - 1);
				nested_loop<D - M>(start, end, [&](const std::array<size_t, D - M> &offset) {
					LocalData<M> this_ghost = data.getGhostSliceOn(f, offset);
					nested_loop<M>(this_ghost.getStart(), this_ghost.getEnd(), [&](const std::array<int, M> &coord) { this_ghost[coord] = 0; });
				});
			}
		}
	}

	void zeroGhostCells(std::shared_ptr<const Vector<D>> u) const
	{
		for (const PatchInfo<D> &pinfo : domain->getPatchInfoVector()) {
			for (auto &this_patch : u->getLocalDatas(pinfo.local_index)) {
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
				if constexpr (D >= 2) {
					ghost_local_data_infos.template get<0>() = GhostLocalDataInfo<0>(*domain);
					enumerateCalls<0>(local_calls.template get<0>(), rank_to_remote_call_sets);
				}
			case GhostFillingType::Edges:
				if constexpr (D == 3) {
					ghost_local_data_infos.template get<1>() = GhostLocalDataInfo<1>(*domain);
					enumerateCalls<1>(local_calls.template get<1>(), rank_to_remote_call_sets);
				}
			case GhostFillingType::Faces:
				ghost_local_data_infos.template get<D - 1>() = GhostLocalDataInfo<D - 1>(*domain);
				enumerateCalls<D - 1>(local_calls.template get<D - 1>(), rank_to_remote_call_sets);
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
		zeroGhostCells(u);

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
		processLocalFills<D - 1>(u);

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
