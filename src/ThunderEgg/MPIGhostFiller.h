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

#include <ThunderEgg/GhostFiller.h>
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
template <size_t D> class MPIGhostFiller : public GhostFiller<D>
{
	private:
	/**
	 * @brief last three ints are the local index for the local patch, and index in the buffer
	 */
	using RemoteCall = std::tuple<std::shared_ptr<const PatchInfo<D>>, const Side<D>, const NbrType,
	                              const Orthant<D>, int, size_t>;
	/**
	 * @brief last two ints are local indexes for the local and nbr patch
	 */
	using LocalCall = std::tuple<std::shared_ptr<const PatchInfo<D>>, const Side<D>, const NbrType,
	                             const Orthant<D>, int, int>;

	/**
	 * @brief A vector of RemoteCall deques, one sperate deque for each rank
	 */
	std::vector<std::deque<RemoteCall>> remote_calls;
	/**
	 * @brief deque of local calls to be made
	 */
	std::deque<LocalCall> local_calls;
	/**
	 * @brief vector of deques for incoming ghost cells, one deque for each rank
	 *
	 * the deques contain a tuple with the local index of the patch, the side that the ghost cells
	 * are on, and the offset in the buffer for those ghost cells
	 */
	std::vector<std::deque<std::tuple<int, Side<D>, size_t>>> incoming_ghosts;
	/**
	 * @brief vectors ranks, the position of the ranks correlate with other vectors.
	 */
	std::vector<size_t> index_rank_map;
	/**
	 * @brief lengths of send buffers for each rank
	 */
	std::vector<size_t> send_buff_lengths;
	/**
	 * @brief lengths of recv buffers for each rank
	 */
	std::vector<size_t> recv_buff_lengths;

	/**
	 * @brief Fill the ghost cells for the neighboring patch
	 *
	 * @param pinfo the patch that ghost cells are being filled from
	 * @param local_data the local data for patch that ghost cells are being filled from
	 * @param nbr_data  the local data for the neighboring patch, where ghost cells are being
	 * filled.
	 * @param side the sided that the neighboring patch is on
	 * @param nbr_type the type of neighbor
	 * @param orthant the orthant that the neighbors ghost cells lie on
	 */
	virtual void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                       const LocalData<D>                  local_data,
	                                       const LocalData<D> nbr_data, const Side<D> side,
	                                       const NbrType    nbr_type,
	                                       const Orthant<D> orthant) const = 0;

	/**
	 * @brief Perform any on this patches ghost cells.
	 *
	 * This may be necessary on some schemes because it needs data from the patch itself, not just
	 * the neighboring patch
	 *
	 * @param pinfo the patch
	 * @param local_data the LocalData for the patch
	 */
	virtual void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                         const LocalData<D> local_data) const = 0;

	/**
	 * @brief Get the LocalData object for the buffer
	 *
	 * @param buffer_ptr pointer to the ghost cells position in the buffer
	 * @param side  the side that the ghost cells are on
	 * @return LocalData<D> the LocalData object
	 */
	LocalData<D> getLocalDataForBuffer(double *buffer_ptr, const Side<D> side) const
	{
		auto ns              = domain->getNs();
		int  num_ghost_cells = domain->getNumGhostCells();
		// determine striding
		std::array<int, D> strides;
		strides[0] = 1;
		for (size_t i = 1; i < D; i++) {
			if (i == side.getAxisIndex() + 1) {
				strides[i] = num_ghost_cells * strides[i - 1];
			} else {
				strides[i] = ns[i - 1] * strides[i - 1];
			}
		}
		// transform buffer ptr so that it points to first non-ghost cell
		double *transformed_buffer_ptr;
		if (side.isLowerOnAxis()) {
			transformed_buffer_ptr = buffer_ptr - (-num_ghost_cells) * strides[side.getAxisIndex()];
		} else {
			transformed_buffer_ptr
			= buffer_ptr - ns[side.getAxisIndex()] * strides[side.getAxisIndex()];
		}

		LocalData<D> buffer_data(transformed_buffer_ptr, strides, ns, num_ghost_cells);
		return buffer_data;
	}

	protected:
	/**
	 * @brief The domain that this ghostfiller operates on
	 */
	std::shared_ptr<const Domain<D>> domain;
	/**
	 * @brief Number of sides to address
	 */
	int side_cases;

	public:
	/**
	 * @brief Construct a new MPIGhostFiller object
	 *
	 * @param domain_in  the domain being used
	 * @param side_cases_in  the number of side cases to address
	 */
	MPIGhostFiller(std::shared_ptr<const Domain<D>> domain_in, int side_cases_in)
	: domain(domain_in), side_cases(side_cases_in)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		// rank nbr_id side orthant nbrtype pinfo local_index
		std::set<std::tuple<int, int, const Side<D>, const Orthant<D>, const NbrType,
		                    std::shared_ptr<const PatchInfo<D>>, int>>
		remote_call_set;
		// rank id side local index
		std::set<std::tuple<int, int, const Side<D>, int>> incoming_ghost_set;

		std::array<size_t, D> axis_ghost_lengths;
		for (size_t axis = 0; axis < D; axis++) {
			size_t length = 1;
			for (size_t i = 0; i < D; i++) {
				if (i == axis) {
					length *= domain->getNumGhostCells();
				} else {
					length *= domain->getNs()[i];
				}
			}
			axis_ghost_lengths[axis] = length;
		}
		std::set<int> ranks;
		for (auto pinfo : domain->getPatchInfoVector()) {
			for (Side<D> s : Side<D>::getValues()) {
				if (pinfo->hasNbr(s)) {
					switch (pinfo->getNbrType(s)) {
						case NbrType::Normal: {
							auto nbrinfo = pinfo->getNormalNbrInfo(s);
							if (nbrinfo.rank == rank) {
								local_calls.emplace_back(pinfo, s, NbrType::Normal,
								                         Orthant<D>::null(), pinfo->local_index,
								                         nbrinfo.local_index);
							} else {
								ranks.insert(nbrinfo.rank);
								remote_call_set.emplace(nbrinfo.rank, nbrinfo.id, s.opposite(),
								                        Orthant<D>::null(), NbrType::Normal, pinfo,
								                        pinfo->local_index);
								incoming_ghost_set.emplace(nbrinfo.rank, pinfo->id, s,
								                           pinfo->local_index);
							}
						} break;
						case NbrType::Fine: {
							auto nbrinfo  = pinfo->getFineNbrInfo(s);
							auto orthants = Orthant<D>::getValuesOnSide(s);
							for (size_t i = 0; i < orthants.size(); i++) {
								if (nbrinfo.ranks[i] == rank) {
									local_calls.emplace_back(pinfo, s, NbrType::Fine, orthants[i],
									                         pinfo->local_index,
									                         nbrinfo.local_indexes[i]);
								} else {
									ranks.insert(nbrinfo.ranks[i]);
									remote_call_set.emplace(
									nbrinfo.ranks[i], nbrinfo.ids[i], s.opposite(), orthants[i],
									NbrType::Fine, pinfo, pinfo->local_index);
									incoming_ghost_set.emplace(nbrinfo.ranks[i], pinfo->id, s,
									                           pinfo->local_index);
								}
							}
						} break;
						case NbrType::Coarse: {
							auto nbrinfo = pinfo->getCoarseNbrInfo(s);
							auto orthant = Orthant<D>::getValuesOnSide(
							s.opposite())[nbrinfo.orth_on_coarse.getIndex()];
							if (nbrinfo.rank == rank) {
								local_calls.emplace_back(pinfo, s, NbrType::Coarse, orthant,
								                         pinfo->local_index, nbrinfo.local_index);
							} else {
								ranks.insert(nbrinfo.rank);
								remote_call_set.emplace(nbrinfo.rank, nbrinfo.id, s.opposite(),
								                        orthant, NbrType::Coarse, pinfo,
								                        pinfo->local_index);
								incoming_ghost_set.emplace(nbrinfo.rank, pinfo->id, s,
								                           pinfo->local_index);
							}
						} break;
					}
				}
			}
		}
		std::map<int, size_t> rank_index_map;
		index_rank_map.reserve(ranks.size());
		int curr_index = 0;
		for (int rank : ranks) {
			rank_index_map[rank] = curr_index;
			index_rank_map.push_back(rank);
			curr_index++;
		}
		std::tuple<int, Side<D>> prev_id_side;
		send_buff_lengths.resize(ranks.size());
		remote_calls.resize(ranks.size());
		for (auto call : remote_call_set) {
			int  local_buffer_index = rank_index_map[std::get<0>(call)];
			int  id                 = std::get<1>(call);
			auto side               = std::get<2>(call).opposite();
			auto orthant            = std::get<3>(call);
			auto nbr_type           = std::get<4>(call);
			auto pinfo              = std::get<5>(call);
			auto local_index        = std::get<6>(call);

			size_t offset = send_buff_lengths[local_buffer_index];

			// calculate length in buffer need for ghost cells
			size_t length = axis_ghost_lengths[side.getAxisIndex()];
			// add length to buffer length
			// if its the same side of the patch, resuse the previous buffer space
			if (std::make_tuple(id, side) == prev_id_side) {
				offset -= length;
			} else {
				send_buff_lengths[local_buffer_index] += length;
			}

			remote_calls[local_buffer_index].emplace_back(pinfo, side, nbr_type, orthant,
			                                              local_index, offset);
			prev_id_side = std::make_tuple(id, side);
		}
		recv_buff_lengths.resize(ranks.size());
		incoming_ghosts.resize(ranks.size());
		for (auto t : incoming_ghost_set) {
			int     local_buffer_index = rank_index_map[std::get<0>(t)];
			Side<D> side               = std::get<2>(t);
			int     local_index        = std::get<3>(t);

			// add length for ghosts to buffer length
			size_t length = axis_ghost_lengths[side.getAxisIndex()];
			size_t offset = recv_buff_lengths[local_buffer_index];
			recv_buff_lengths[local_buffer_index] += length;

			// add ghost to incoming ghosts
			incoming_ghosts[local_buffer_index].emplace_back(local_index, side, offset);
		}
	}
	/**
	 * @brief Fill ghost cells on a vector
	 *
	 * @param u  the vector
	 */
	void fillGhost(std::shared_ptr<const Vector<D>> u) const
	{
		// zero out ghost cells
		for (auto pinfo : domain->getPatchInfoVector()) {
			const LocalData<2> this_patch = u->getLocalData(pinfo->local_index);
			for (Side<2> s : Side<2>::getValues()) {
				if (pinfo->hasNbr(s)) {
					for (int i = 0; i < pinfo->num_ghost_cells; i++) {
						auto this_ghost = this_patch.getGhostSliceOnSide(s, i + 1);
						nested_loop<1>(
						this_ghost.getStart(), this_ghost.getEnd(),
						[&](const std::array<int, 1> &coord) { this_ghost[coord] = 0; });
					}
				}
			}
		}

		// allocate recv buffers and post recvs
		std::vector<std::vector<double>> recv_buffers(recv_buff_lengths.size());
		for (size_t i = 0; i < recv_buff_lengths.size(); i++) {
			recv_buffers[i].resize(recv_buff_lengths[i]);
		}
		std::vector<MPI_Request> recv_requests(recv_buff_lengths.size());
		for (size_t i = 0; i < recv_requests.size(); i++) {
			MPI_Irecv(recv_buffers[i].data(), recv_buffers[i].size(), MPI_DOUBLE, index_rank_map[i],
			          0, MPI_COMM_WORLD, &recv_requests[i]);
		}

		// allocate send buffers
		std::vector<std::vector<double>> out_buffers(send_buff_lengths.size());
		for (size_t i = 0; i < send_buff_lengths.size(); i++) {
			out_buffers[i].resize(send_buff_lengths[i]);
		}
		// perform fill operations for mpi buffers
		std::vector<MPI_Request> send_requests(send_buff_lengths.size());
		for (size_t i = 0; i < remote_calls.size(); i++) {
			for (const RemoteCall &call : remote_calls[i]) {
				auto    pinfo         = std::get<0>(call);
				auto    side          = std::get<1>(call);
				auto    nbr_type      = std::get<2>(call);
				auto    orthant       = std::get<3>(call);
				auto    local_data    = u->getLocalData(std::get<4>(call));
				size_t  buffer_offset = std::get<5>(call);
				double *buffer_ptr    = out_buffers[i].data() + buffer_offset;

				// create a LocalData object for the buffer
				LocalData<D> buffer_data = getLocalDataForBuffer(buffer_ptr, side.opposite());

				// make the call
				fillGhostCellsForNbrPatch(pinfo, local_data, buffer_data, side, nbr_type, orthant);
			}
			MPI_Isend(out_buffers[i].data(), out_buffers[i].size(), MPI_DOUBLE, index_rank_map[i],
			          0, MPI_COMM_WORLD, &send_requests[i]);
		}
		// perform local operations
		for (auto pinfo : domain->getPatchInfoVector()) {
			auto data = u->getLocalData(pinfo->local_index);
			fillGhostCellsForLocalPatch(pinfo, data);
		}
		for (const LocalCall &call : local_calls) {
			auto pinfo      = std::get<0>(call);
			auto side       = std::get<1>(call);
			auto nbr_type   = std::get<2>(call);
			auto orthant    = std::get<3>(call);
			auto local_data = u->getLocalData(std::get<4>(call));
			auto nbr_data   = u->getLocalData(std::get<5>(call));
			fillGhostCellsForNbrPatch(pinfo, local_data, nbr_data, side, nbr_type, orthant);
		}
		// add in remote values as they come in
		for (size_t i = 0; i < recv_requests.size(); i++) {
			int finished_index;
			MPI_Waitany(recv_requests.size(), recv_requests.data(), &finished_index,
			            MPI_STATUS_IGNORE);
			for (auto t : incoming_ghosts[finished_index]) {
				int     local_index   = std::get<0>(t);
				Side<D> side          = std::get<1>(t);
				size_t  buffer_offset = std::get<2>(t);

				const LocalData<D> local_data = u->getLocalData(local_index);
				double *           buffer_ptr = recv_buffers[finished_index].data() + buffer_offset;
				LocalData<D>       buffer_data = getLocalDataForBuffer(buffer_ptr, side);
				for (int ig = 0; ig < domain->getNumGhostCells(); ig++) {
					LocalData<D - 1> local_slice  = local_data.getGhostSliceOnSide(side, ig + 1);
					LocalData<D - 1> buffer_slice = buffer_data.getGhostSliceOnSide(side, ig + 1);
					nested_loop<D - 1>(local_slice.getStart(), local_slice.getEnd(),
					                   [&](const std::array<int, D - 1> &coord) {
						                   local_slice[coord] += buffer_slice[coord];
					                   });
				}
			}
		}

		// wait for sends for finish
		MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUS_IGNORE);
	}
};
} // namespace ThunderEgg
#endif