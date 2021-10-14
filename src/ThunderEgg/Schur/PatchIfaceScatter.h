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

#ifndef THUNDEREGG_SCHUR_PATCHIFACESCATTER_H
#define THUNDEREGG_SCHUR_PATCHIFACESCATTER_H
/**
 * @file
 *
 * @brief PatchIfaceScatter class
 */
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>
#include <ThunderEgg/Vector.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief Scatters between a global Schur compliment vector and a local patch iface vector
 *
 * The scatters functions are split with a Start and Finish, this allows for local computation to
 * occur while the communicating
 *
 * @tparam D the number of cartesian dimensions on a patch
 */
template <int D> class PatchIfaceScatter
{
	private:
	class StatePrivate
	{
		public:
		std::vector<std::vector<double>> send_buffers;
		std::vector<std::vector<double>> recv_buffers;
		std::vector<MPI_Request>         send_requests;
		std::vector<MPI_Request>         recv_requests;
		bool                             communicating = true;
		/**
		 * @brief The global vector passed to scatterStart
		 */
		const Vector<D - 1> *curr_global_vector = nullptr;
		/**
		 * @brief The local vector passed to scatterStart
		 */
		const Vector<D - 1> *curr_local_vector = nullptr;

		StatePrivate(size_t send_buffers_size, size_t recv_buffers_size)
		: send_buffers(send_buffers_size),
		  recv_buffers(recv_buffers_size),
		  send_requests(send_buffers_size),
		  recv_requests(recv_buffers_size)
		{
		}
		~StatePrivate()
		{
			if (communicating) {
				MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
				MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
			}
		}
	};

	class State
	{
		public:
		std::shared_ptr<StatePrivate> ptr;
		State(size_t num_send, size_t num_recv) : ptr(new StatePrivate(num_send, num_recv)) {}
	};

	public:
	/**
	 * @brief the number of cells in each direction of the interface
	 */
	std::array<int, D - 1> lengths;
	/**
	 * @brief the number of interfaces in the local vector
	 */
	int num_local_patch_ifaces;
	/**
	 * @brief the number of interfaces in the local vector
	 */
	int iface_stride;
	/**
	 * @brief number of MPI sends
	 */
	int num_sends;
	/**
	 * @brief MPI Send ranks
	 */
	std::vector<int> send_ranks;
	/**
	 * @brief the local indexes of the local vector to send
	 */
	std::vector<std::vector<int>> send_local_indexes;
	/**
	 * @brief number of MPI recvs
	 */
	int num_recvs;
	/**
	 * @brief MPI recv ranks
	 */
	std::vector<int> recv_ranks;
	/**
	 * @brief local indexes of the local vector to receive
	 */
	std::vector<std::vector<int>> recv_local_indexes;

	/**
	 * @brief Set the incoming buffer maps (incoming_rank_to_local_indexes) and set
	 * num_local_patch_ifaces
	 *
	 * @param iface_domain the InterfaceDomain
	 */
	void setIncomingBufferMapsAndDetermineLocalVectorSize(const InterfaceDomain<D> &iface_domain)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		std::map<int, std::set<std::pair<int, int>>> incoming_ranks_to_id_local_index_pairs;

		for (auto piinfo : iface_domain.getPatchIfaceInfos()) {
			for (Side<D> s : Side<D>::getValues()) {
				if (piinfo->pinfo.hasNbr(s)) {
					auto iface_info = piinfo->getIfaceInfo(s);
					if (iface_info->rank != rank) {
						incoming_ranks_to_id_local_index_pairs[iface_info->rank].emplace(iface_info->id, iface_info->patch_local_index);
					}
				}
			}
		}
		int local_vector_size = 0;
		for (auto rank_to_id_local_index_pairs : incoming_ranks_to_id_local_index_pairs) {
			local_vector_size += rank_to_id_local_index_pairs.second.size();
		}
		local_vector_size += iface_domain.getNumLocalInterfaces();
		num_local_patch_ifaces = local_vector_size;

		num_recvs = incoming_ranks_to_id_local_index_pairs.size();
		recv_ranks.resize(num_recvs);
		recv_local_indexes.resize(num_recvs);

		int recv_index = 0;
		for (const auto &rank_to_id_local_index_pairs : incoming_ranks_to_id_local_index_pairs) {
			recv_ranks[recv_index]               = rank_to_id_local_index_pairs.first;
			std::vector<int> &local_index_vector = recv_local_indexes[recv_index];

			local_index_vector.reserve(rank_to_id_local_index_pairs.second.size());

			for (auto id_local_index_pair : rank_to_id_local_index_pairs.second) {
				local_index_vector.push_back(id_local_index_pair.second);
			}
			recv_index++;
		}
	}
	/**
	 * @brief Set the outgoing buffer maps (outgoing_rank_to_local_indexes)
	 *
	 * @param iface_domain the InterfaceDomain
	 */
	void setOutgoingBufferMaps(const InterfaceDomain<D> &iface_domain)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		std::map<int, std::set<std::pair<int, int>>> outgoing_ranks_to_id_local_index_pairs;
		for (auto iface : iface_domain.getInterfaces()) {
			for (auto patch : iface->patches) {
				if ((patch.type.isNormal() || patch.type.isFineToFine() || patch.type.isCoarseToCoarse()) && patch.piinfo->pinfo.rank != rank) {
					outgoing_ranks_to_id_local_index_pairs[patch.piinfo->pinfo.rank].emplace(iface->id, iface->local_index);
				}
			}
		}

		num_sends = outgoing_ranks_to_id_local_index_pairs.size();
		send_ranks.resize(num_sends);
		send_local_indexes.resize(num_sends);

		int send_index = 0;
		for (const auto &rank_to_id_local_index_pairs : outgoing_ranks_to_id_local_index_pairs) {
			send_ranks[send_index]               = rank_to_id_local_index_pairs.first;
			std::vector<int> &local_index_vector = send_local_indexes[send_index];

			local_index_vector.reserve(rank_to_id_local_index_pairs.second.size());

			for (auto id_local_index_pair : rank_to_id_local_index_pairs.second) {
				local_index_vector.push_back(id_local_index_pair.second);
			}
			send_index++;
		}
	}
	/**
	 * @brief Initialize the mpi buffers
	 */
	State initializeMPIBuffers() const
	{
		State state(send_ranks.size(), recv_ranks.size());

		for (int send_index = 0; send_index < num_sends; send_index++) {
			state.ptr->send_buffers[send_index].resize(send_local_indexes[send_index].size() * iface_stride);
		}

		for (int recv_index = 0; recv_index < num_recvs; recv_index++) {
			state.ptr->recv_buffers[recv_index].resize(recv_local_indexes[recv_index].size() * iface_stride);
		}

		return state;
	}

	public:
	/**
	 * @brief Construct a new PatchIfaceScatter object
	 *
	 * @param iface_domain the InterfaceDomain
	 */
	explicit PatchIfaceScatter(const InterfaceDomain<D> &iface_domain)
	{
		std::array<int, D> ns = iface_domain.getDomain().getNs();
		for (int i = 1; i < D; i++) {
			if (ns[0] != ns[i]) {
				throw RuntimeError("Cannot form Schur compliment vector for Domain with non-square patches");
			}
		}

		lengths.fill(ns[0]);
		iface_stride = std::pow(ns[0], D - 1);
		setIncomingBufferMapsAndDetermineLocalVectorSize(iface_domain);
		setOutgoingBufferMaps(iface_domain);
	}
	/**
	 * @brief Get a nw local patch iface vector
	 *
	 * @return std::shared_ptr<Vector<D - 1>> the new vector
	 */
	std::shared_ptr<Vector<D - 1>> getNewLocalPatchIfaceVector() const
	{
		return std::make_shared<Vector<D - 1>>(Communicator(MPI_COMM_SELF), lengths, 1, num_local_patch_ifaces, 0);
	}
	/**
	 * @brief Start the scatter from the global Schur compliment vector to the local patch iface
	 * vector
	 *
	 * Interfaces that are local to this processor will be copied to the local_patch_iface_vector at
	 * the end of this call
	 *
	 * Will throw an exception if any communcation is in progress
	 *
	 * @param global_vector the global Schur compliment vector
	 * @param local_patch_iface_vector the the local patch iface vector
	 */
	State scatterStart(const Vector<D - 1> &global_vector, Vector<D - 1> &local_patch_iface_vector) const
	{
		for (int i = 0; i < global_vector.getNumLocalPatches(); i++) {
			auto global_data = global_vector.getComponentView(0, i);
			auto local_data  = local_patch_iface_vector.getComponentView(0, i);
			Loop::Nested<D - 1>(
			local_data.getStart(), local_data.getEnd(), [&](const std::array<int, D - 1> &coord) { local_data[coord] = global_data[coord]; });
		}

		State state = initializeMPIBuffers();

		for (int recv_index = 0; recv_index < num_recvs; recv_index++) {
			MPI_Irecv(state.ptr->recv_buffers[recv_index].data(),
			          state.ptr->recv_buffers[recv_index].size(),
			          MPI_DOUBLE,
			          recv_ranks[recv_index],
			          0,
			          MPI_COMM_WORLD,
			          &state.ptr->recv_requests[recv_index]);
		}

		for (int send_index = 0; send_index < num_sends; send_index++) {
			std::vector<double> &buffer = state.ptr->send_buffers[send_index];

			int buffer_index = 0;
			for (int local_index : send_local_indexes[send_index]) {
				auto local_data = global_vector.getComponentView(0, local_index);
				Loop::OverInteriorIndexes<D - 1>(local_data, [&](const std::array<int, D - 1> &coord) {
					buffer[buffer_index] = local_data[coord];
					buffer_index++;
				});
			}

			MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, send_ranks[send_index], 0, MPI_COMM_WORLD, &state.ptr->send_requests[send_index]);
		}

		for (int local_iface = 0; local_iface < global_vector.getNumLocalPatches(); local_iface++) {
			auto global_data = global_vector.getComponentView(0, local_iface);
			auto local_data  = local_patch_iface_vector.getComponentView(0, local_iface);
			Loop::OverInteriorIndexes<D - 1>(local_data, [&](const std::array<int, D - 1> &coord) { local_data[coord] = global_data[coord]; });
		}

		state.ptr->curr_global_vector = &global_vector;
		state.ptr->curr_local_vector  = &local_patch_iface_vector;
		return state;
	}
	/**
	 * @brief Finish the scatter from the global Schur compliment vector to the local patch iface
	 * vector
	 *
	 * Will throw an exception if different vectors are passed
	 * from when scatterStart was called
	 *
	 * @param global_vector the global Schur compliment vector
	 * @param local_patch_iface_vector the the local patch iface vector
	 */
	void scatterFinish(const State &state, const Vector<D - 1> &global_vector, Vector<D - 1> &local_patch_iface_vector) const
	{
		if (&global_vector != state.ptr->curr_global_vector || &local_patch_iface_vector != state.ptr->curr_local_vector) {
			throw RuntimeError("Different vectors were passed ot scatterFinish than were passed to scatterStart");
		}

		for (int i = 0; i < num_recvs; i++) {
			MPI_Status status;
			int        recv_index;
			MPI_Waitany(num_recvs, state.ptr->recv_requests.data(), &recv_index, &status);

			std::vector<double> &buffer = state.ptr->recv_buffers[recv_index];

			int buffer_index = 0;
			for (int local_index : recv_local_indexes[recv_index]) {
				auto local_data = local_patch_iface_vector.getComponentView(0, local_index);
				Loop::OverInteriorIndexes<D - 1>(local_data, [&](const std::array<int, D - 1> &coord) {
					local_data[coord] = buffer[buffer_index];
					buffer_index++;
				});
			}
		}

		MPI_Waitall(num_sends, state.ptr->send_requests.data(), MPI_STATUSES_IGNORE);

		state.ptr->communicating = false;
	}
};
extern template class PatchIfaceScatter<2>;
extern template class PatchIfaceScatter<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
