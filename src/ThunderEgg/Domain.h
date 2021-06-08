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

#ifndef THUNDEREGG_DOMAIN_H
#define THUNDEREGG_DOMAIN_H
#include <ThunderEgg/PatchInfo.h>
#include <ThunderEgg/Timer.h>
#include <ThunderEgg/Vector.h>
#include <cmath>
#include <deque>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace ThunderEgg
{
/**
 * @brief Uses a collection of PatchInfo objects to represent the domain of the problem.
 *
 * Each patch passed to the constructor needs to have the following complete information:
 * 	* Each patch within a domain have to have a unique id set.
 *  * Each patch needs to have the neighboring patch information filled out, with the id of neighbor
 * patches.
 *
 *  When passed to the constructor, the constructor will create an indexing for the
 * This class mainly manages a set of patches that makes up the domain. It is responsible for
 * setting up the indexing of the domains, which is used in the rest of the ThunderEgg library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class Domain
{
	private:
	/**
	 * @brief The id of the domain
	 */
	int id = -1;
	/**
	 * @brief The number of cells in each direction
	 */
	std::array<int, D> ns;
	/**
	 * @brief number of ghost cells on each side of the patch
	 */
	int num_ghost_cells;
	/**
	 * @brief The number of cells in a patch
	 */
	int num_cells_in_patch;
	/**
	 * @brief The number of cells(including ghost cells) in a patch
	 */
	int num_cells_in_patch_with_ghost;
	/**
	 * @brief Vector of PatchInfo pointers where index in the vector corresponds to the patch's
	 * local index
	 */
	std::vector<PatchInfo<D>> pinfos;
	/**
	 * @brief The global number of patches
	 */
	int global_num_patches = 1;
	/**
	 * @brief The timer
	 */
	mutable std::shared_ptr<Timer> timer;

	/**
	 * @brief Give the patches local indexes.
	 */
	void indexPatchesLocal()
	{
		// index patches
		int                curr_index = 0;
		std::map<int, int> id_to_local_index;
		for (auto &pinfo : pinfos) {
			pinfo.local_index           = curr_index;
			id_to_local_index[pinfo.id] = pinfo.local_index;
			curr_index++;
		}

		// set local index in nbrinfo objects
		for (auto &pinfo : pinfos) {
			pinfo.setNeighborLocalIndexes(id_to_local_index);
		}
	}
	/**
	 * @brief Give the patches global indexes
	 */
	void indexPatchesGlobal()
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// get starting global index
		int num_local_patches = (int) pinfos.size();
		int curr_global_index;
		MPI_Scan(&num_local_patches, &curr_global_index, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		curr_global_index -= num_local_patches;

		// index the patches
		std::map<int, int> id_to_global_index;
		for (auto &pinfo : pinfos) {
			pinfo.global_index           = curr_global_index;
			id_to_global_index[pinfo.id] = pinfo.global_index;
			curr_global_index++;
		}

		std::map<int, std::set<std::pair<int, int>>> ranks_to_ids_and_global_indexes_outgoing;
		std::map<int, std::set<int>>                 ranks_to_ids_incoming;
		for (auto &pinfo : pinfos) {
			auto ranks = pinfo.getNbrRanks();
			auto ids   = pinfo.getNbrIds();
			for (size_t idx = 0; idx < ranks.size(); idx++) {
				int nbr_id   = ids[idx];
				int nbr_rank = ranks[idx];
				if (nbr_rank != rank) {
					ranks_to_ids_and_global_indexes_outgoing[nbr_rank].insert(std::make_pair(pinfo.id, pinfo.global_index));
					ranks_to_ids_incoming[nbr_rank].insert(nbr_id);
				}
			}
		}

		// prepare to recieve data
		std::vector<MPI_Request> recv_requests;

		// allocate incoming vectors and post recvs
		std::map<int, std::vector<int>> rank_to_incoming_data;
		for (const auto &pair : ranks_to_ids_incoming) {
			int               source_rank   = pair.first;
			std::vector<int> &incoming_data = rank_to_incoming_data[source_rank];
			incoming_data.resize(pair.second.size());

			MPI_Request request;
			MPI_Irecv(incoming_data.data(), (int) incoming_data.size(), MPI_INT, source_rank, 0, MPI_COMM_WORLD, &request);
			recv_requests.push_back(request);
		}

		// prepare outgoing vector of data and send it
		std::vector<MPI_Request> send_requests;

		// post sends
		std::map<int, std::vector<int>> rank_to_outgoing_data;
		for (const auto &pair : ranks_to_ids_and_global_indexes_outgoing) {
			int dest_rank = pair.first;
			// allocate and fill vector
			std::vector<int> &data = rank_to_outgoing_data[dest_rank];
			data.reserve(pair.second.size());
			for (const auto &id_and_global_index : pair.second) {
				data.push_back(id_and_global_index.second);
			}
			MPI_Request request;
			MPI_Isend(data.data(), (int) data.size(), MPI_INT, dest_rank, 0, MPI_COMM_WORLD, &request);
			send_requests.push_back(request);
		}

		// add global indexes to map as recvs come in
		for (size_t i = 0; i < recv_requests.size(); i++) {
			MPI_Status status;
			int        request_index;
			MPI_Waitany((int) recv_requests.size(), recv_requests.data(), &request_index, &status);

			int                     source_rank  = status.MPI_SOURCE;
			const std::set<int> &   incoming_ids = ranks_to_ids_incoming[source_rank];
			const std::vector<int> &data         = rank_to_incoming_data[source_rank];

			auto curr_id                = incoming_ids.cbegin();
			auto curr_global_index_iter = data.cbegin();
			while (curr_id != incoming_ids.cend()) {
				id_to_global_index[*curr_id] = *curr_global_index_iter;
				curr_id++;
				curr_global_index_iter++;
			}
		}

		// wait for all the sends to finsh
		MPI_Waitall((int) send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

		// update global indexes in nbrinfo objects
		for (auto &pinfo : pinfos) {
			pinfo.setNeighborGlobalIndexes(id_to_global_index);
		}
	}

	public:
	/**
	 * @brief Construct a new Domain object
	 *
	 * @tparam InputIterator the iterator for PatchInfo objects
	 * @param id the id of the domain should be unique within a multigrid cycle
	 * @param ns the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @param first_pinfo start iterator for PatchInfo objects
	 * @param last_pinfo end iterator for PatchInfo objects
	 */
	template <class InputIterator>
	Domain(int id, std::array<int, D> ns, int num_ghost_cells, InputIterator first_pinfo, InputIterator last_pinfo)
	: id(id),
	  ns(ns),
	  num_ghost_cells(num_ghost_cells),
	  pinfos(first_pinfo, last_pinfo)
	{
		num_cells_in_patch            = 1;
		num_cells_in_patch_with_ghost = 1;
		for (size_t i = 0; i < D; i++) {
			num_cells_in_patch *= ns[i];
			num_cells_in_patch_with_ghost *= (ns[i] + 2 * num_ghost_cells);
		}

		int num_local_domains = pinfos.size();
		MPI_Allreduce(&num_local_domains, &global_num_patches, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		indexPatchesLocal();
		indexPatchesGlobal();
	}
	/**
	 * @brief Get a vector of PatchInfo pointers where index in the vector corresponds to the
	 * patch's local index
	 */
	const std::vector<PatchInfo<D>> &getPatchInfoVector() const
	{
		return pinfos;
	}
	/**
	 * @brief Get the number of cells in each direction
	 *
	 */
	const std::array<int, D> &getNs() const
	{
		return ns;
	}
	/**
	 * @brief Get the number of global patches
	 */
	int getNumGlobalPatches() const
	{
		return global_num_patches;
	}
	/**
	 * @brief Get the number of local patches
	 */
	int getNumLocalPatches() const
	{
		return (int) pinfos.size();
	}
	/**
	 * @brief get the number of global cells
	 */
	int getNumGlobalCells() const
	{
		return global_num_patches * num_cells_in_patch;
	}
	/**
	 * @brief Get get the number of local cells
	 */
	int getNumLocalCells() const
	{
		return ((int) pinfos.size()) * num_cells_in_patch;
	}
	/**
	 * @brief Get get the number of local cells (including ghost cells)
	 */
	int getNumLocalCellsWithGhost() const
	{
		return ((int) pinfos.size()) * num_cells_in_patch_with_ghost;
	}
	/**
	 * @brief Get the number of cells in a patch
	 */
	int getNumCellsInPatch() const
	{
		return num_cells_in_patch;
	}
	/**
	 * @brief get the number of ghost cell on each side of a patch
	 */
	int getNumGhostCells() const
	{
		return num_ghost_cells;
	}
	/**
	 * @brief Get the volume of the domain.
	 *
	 * For 2D, this will be the area.
	 */
	double volume() const
	{
		double sum = 0;
		for (auto &pinfo : pinfos) {
			double patch_vol = 1;
			for (size_t i = 0; i < D; i++) {
				patch_vol *= pinfo.spacings[i] * pinfo.ns[i];
			}
			sum += patch_vol;
		}
		double retval;
		MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return retval;
	}
	/**
	 * @brief Integrate a vector over the domain.
	 *
	 * @param u the vector
	 * @return double the result of the integral
	 */
	double integrate(std::shared_ptr<const Vector<D>> u) const
	{
		double sum = 0;

		for (const auto &pinfo : pinfos) {
			for (int c = 0; c < u->getNumComponents(); c++) {
				ComponentView<const double, D> u_data = u->getComponentView(c, pinfo.local_index);

				double patch_sum = 0;
				nested_loop<D>(u_data.getStart(), u_data.getEnd(), [&](std::array<int, D> coord) { patch_sum += u_data[coord]; });

				for (size_t i = 0; i < D; i++) {
					patch_sum *= pinfo.spacings[i];
				}
				sum += patch_sum;
			}
		}
		double retval;
		MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return retval;
	}
	/**
	 * @brief Set the Timer object
	 *
	 * @param timer the timer
	 */
	void setTimer(std::shared_ptr<Timer> timer) const
	{
		this->timer = timer;
		timer->addDomain(id, *this);
	}
	/**
	 * @brief Get the Timer object
	 *
	 * @return std::shared_ptr<Timer> the timer
	 */
	std::shared_ptr<Timer> getTimer() const
	{
		return timer;
	}
	/**
	 * @brief Check if the Domain has a timer associated with it
	 * @return true if the Domain has a timer associated with it
	 */
	bool hasTimer() const
	{
		return timer != nullptr;
	}
	/**
	 * @brief Get the domain's id
	 *
	 * @return int the id
	 */
	int getId() const
	{
		return id;
	}
};

template <int D> void to_json(nlohmann::json &j, const Domain<D> &domain)
{
	for (auto pinfo : domain.getPatchInfoVector()) {
		j.push_back(pinfo);
	}
}

extern template class Domain<2>;
extern template class Domain<3>;
extern template void to_json<2>(nlohmann::json &j, const Domain<2> &domain);
extern template void to_json<3>(nlohmann::json &j, const Domain<3> &domain);
} // namespace ThunderEgg
#endif
