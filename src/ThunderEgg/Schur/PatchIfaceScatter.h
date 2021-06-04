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

#ifndef THUNDEREGG_SCHUR_PATCHIFACESCATTER_H
#define THUNDEREGG_SCHUR_PATCHIFACESCATTER_H
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>
#include <ThunderEgg/ValVector.h>
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
	/**
	 * @brief are we communicating?
	 */
	bool communicating = false;
	/**
	 * @brief The global vector passed to scatterStart
	 */
	std::shared_ptr<const Vector<D - 1>> curr_global_vector;
	/**
	 * @brief The local vector passed to scatterStart
	 */
	std::shared_ptr<const Vector<D - 1>> curr_local_vector;
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
	 * @brief MPI Send requests
	 */
	std::vector<MPI_Request> send_requests;
	/**
	 * @brief the local indexes of the local vector to send
	 */
	std::vector<std::vector<int>> send_local_indexes;
	/**
	 * @brief MPI Send buffers
	 */
	std::vector<std::vector<double>> send_buffers;
	/**
	 * @brief number of MPI recvs
	 */
	int num_recvs;
	/**
	 * @brief MPI recv ranks
	 */
	std::vector<int> recv_ranks;
	/**
	 * @brief MPI recv requests
	 */
	std::vector<MPI_Request> recv_requests;
	/**
	 * @brief local indexes of the local vector to receive
	 */
	std::vector<std::vector<int>> recv_local_indexes;
	/**
	 * @brief MPI recv buffers
	 */
	std::vector<std::vector<double>> recv_buffers;

	/**
	 * @brief Set the incoming buffer maps (incoming_rank_to_local_indexes) and set
	 * num_local_patch_ifaces
	 *
	 * @param iface_domain the InterfaceDomain
	 */
	void setIncomingBufferMapsAndDetermineLocalVectorSize(std::shared_ptr<const InterfaceDomain<D>> iface_domain)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		std::map<int, std::set<std::pair<int, int>>> incoming_ranks_to_id_local_index_pairs;

		for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
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
		local_vector_size += iface_domain->getNumLocalInterfaces();
		num_local_patch_ifaces = local_vector_size;

		num_recvs = incoming_ranks_to_id_local_index_pairs.size();
		recv_ranks.resize(num_recvs);
		recv_requests.resize(num_recvs);
		recv_local_indexes.resize(num_recvs);
		recv_buffers.resize(num_recvs);

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
	void setOutgoingBufferMaps(std::shared_ptr<const InterfaceDomain<D>> iface_domain)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		std::map<int, std::set<std::pair<int, int>>> outgoing_ranks_to_id_local_index_pairs;
		for (auto iface : iface_domain->getInterfaces()) {
			for (auto patch : iface->patches) {
				if ((patch.type.isNormal() || patch.type.isFineToFine() || patch.type.isCoarseToCoarse()) && patch.piinfo->pinfo.rank != rank) {
					outgoing_ranks_to_id_local_index_pairs[patch.piinfo->pinfo.rank].emplace(iface->id, iface->local_index);
				}
			}
		}

		num_sends = outgoing_ranks_to_id_local_index_pairs.size();
		send_ranks.resize(num_sends);
		send_requests.resize(num_sends);
		send_local_indexes.resize(num_sends);
		send_buffers.resize(num_sends);

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
	void initializeMPIBuffers()
	{
		for (int send_index = 0; send_index < num_sends; send_index++) {
			send_buffers[send_index].resize(send_local_indexes[send_index].size() * iface_stride);
		}

		for (int recv_index = 0; recv_index < num_recvs; recv_index++) {
			recv_buffers[recv_index].resize(recv_local_indexes[recv_index].size() * iface_stride);
		}
	}
	/**
	 * @brief Destroy the mpi buffers
	 */
	void destroyMPIBuffers()
	{
		for (int send_index = 0; send_index < num_sends; send_index++) {
			send_buffers[send_index] = std::vector<double>();
		}
		for (int recv_index = 0; recv_index < num_recvs; recv_index++) {
			recv_buffers[recv_index] = std::vector<double>();
		}
	}

	public:
	/**
	 * @brief Construct a new PatchIfaceScatter object
	 *
	 * @param iface_domain the InterfaceDomain
	 */
	explicit PatchIfaceScatter(std::shared_ptr<const InterfaceDomain<D>> iface_domain)
	{
		std::array<int, D> ns = iface_domain->getDomain()->getNs();
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
	PatchIfaceScatter(const PatchIfaceScatter &) = delete;
	PatchIfaceScatter &operator=(const PatchIfaceScatter &) = delete;
	PatchIfaceScatter(PatchIfaceScatter &&) noexcept        = delete;
	PatchIfaceScatter &operator=(PatchIfaceScatter &&) noexcept = delete;
	/**
	 * @brief Destroy the PatchIfaceScatter object
	 */
	~PatchIfaceScatter()
	{
		if (communicating) {
			MPI_Waitall(num_recvs, recv_requests.data(), MPI_STATUSES_IGNORE);
			MPI_Waitall(num_sends, send_requests.data(), MPI_STATUSES_IGNORE);
		}
	}
	/**
	 * @brief Get a nw local patch iface vector
	 *
	 * @return std::shared_ptr<Vector<D - 1>> the new vector
	 */
	std::shared_ptr<Vector<D - 1>> getNewLocalPatchIfaceVector() const
	{
		return std::make_shared<ValVector<D - 1>>(MPI_COMM_SELF, lengths, 0, 1, num_local_patch_ifaces);
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
	void scatterStart(std::shared_ptr<const Vector<D - 1>> global_vector, std::shared_ptr<Vector<D - 1>> local_patch_iface_vector)
	{
		if (communicating) {
			throw RuntimeError("This PatchIfaceScatter is in the middle of communicating");
		}

		for (int i = 0; i < global_vector->getNumLocalPatches(); i++) {
			auto global_data = global_vector->getComponentView(0, i);
			auto local_data  = local_patch_iface_vector->getComponentView(0, i);
			nested_loop<D - 1>(
			local_data.getStart(), local_data.getEnd(), [&](const std::array<int, D - 1> &coord) { local_data[coord] = global_data[coord]; });
		}

		initializeMPIBuffers();

		for (int recv_index = 0; recv_index < num_recvs; recv_index++) {
			MPI_Irecv(recv_buffers[recv_index].data(),
			          recv_buffers[recv_index].size(),
			          MPI_DOUBLE,
			          recv_ranks[recv_index],
			          0,
			          MPI_COMM_WORLD,
			          &recv_requests[recv_index]);
		}

		for (int send_index = 0; send_index < num_sends; send_index++) {
			std::vector<double> &buffer = send_buffers[send_index];

			int buffer_index = 0;
			for (int local_index : send_local_indexes[send_index]) {
				auto local_data = global_vector->getComponentView(0, local_index);
				nested_loop<D - 1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, D - 1> &coord) {
					buffer[buffer_index] = local_data[coord];
					buffer_index++;
				});
			}

			MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, send_ranks[send_index], 0, MPI_COMM_WORLD, &send_requests[send_index]);
		}

		for (int local_iface = 0; local_iface < global_vector->getNumLocalPatches(); local_iface++) {
			auto global_data = global_vector->getComponentView(0, local_iface);
			auto local_data  = local_patch_iface_vector->getComponentView(0, local_iface);
			nested_loop<D - 1>(
			local_data.getStart(), local_data.getEnd(), [&](const std::array<int, D - 1> &coord) { local_data[coord] = global_data[coord]; });
		}

		curr_global_vector = global_vector;
		curr_local_vector  = local_patch_iface_vector;
		communicating      = true;
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
	void scatterFinish(std::shared_ptr<const Vector<D - 1>> global_vector, std::shared_ptr<Vector<D - 1>> local_patch_iface_vector)
	{
		if (global_vector != curr_global_vector || local_patch_iface_vector != curr_local_vector) {
			throw RuntimeError("Different vectors were passed ot scatterFinish than were passed to scatterStart");
		}

		for (int i = 0; i < num_recvs; i++) {
			MPI_Status status;
			int        recv_index;
			MPI_Waitany(num_recvs, recv_requests.data(), &recv_index, &status);

			std::vector<double> &buffer = recv_buffers[recv_index];

			int buffer_index = 0;
			for (int local_index : recv_local_indexes[recv_index]) {
				auto local_data = local_patch_iface_vector->getComponentView(0, local_index);
				nested_loop<D - 1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, D - 1> &coord) {
					local_data[coord] = buffer[buffer_index];
					buffer_index++;
				});
			}
		}

		MPI_Waitall(num_sends, send_requests.data(), MPI_STATUSES_IGNORE);

		curr_global_vector = nullptr;
		curr_local_vector  = nullptr;
		communicating      = false;
		destroyMPIBuffers();
	}
};
extern template class PatchIfaceScatter<2>;
extern template class PatchIfaceScatter<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
