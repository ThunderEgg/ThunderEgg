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
 * @brief Represents the domain of the problem.
 *
 * This class mainly manages a set of patches that makes up the domain. It is responsible for
 * setting up the indexing of the domains, which is used in the rest of the ThunderEgg library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class Domain
{
	public:
	std::vector<int> patch_id_bc_map_vec;

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
	 * @brief Map that goes form patch's id to the PatchInfo pointer
	 */
	std::map<int, std::shared_ptr<PatchInfo<D>>> pinfo_id_map;
	/**
	 * @brief Vector of PatchInfo pointers where index in the vector corresponds to the patch's
	 * local index
	 */
	std::vector<std::shared_ptr<PatchInfo<D>>> pinfo_vector;
	/**
	 * @brief The global number of patches
	 */
	int global_num_patches = 1;
	/**
	 * @brief local number of boundary conditions in vector
	 */
	int local_num_bc = 4;
	/**
	 * @brief Vector of patch ids. Index corresponds to the patch's local index.
	 */
	std::vector<int> patch_id_map_vec;
	/**
	 * @brief Vector of patch global_indexes. Index corresponds to the patch's global index.
	 */
	std::vector<int> patch_global_index_map_vec;
	/**
	 * @brief vector of neighboring ranks with incoming data in order of how they appear in the
	 * vector
	 */
	std::vector<int> in_off_proc_rank_vec;
	/**
	 * @brief size of incoming data from neighbors
	 */
	std::vector<int> in_off_proc_rank_size_vec;
	/**
	 * @brief map of id to local index in incoming vector
	 */
	std::vector<int> in_off_proc_map_vec;
	/**
	 * @brief vector of neighboring ranks with outgoing data in order of how they appear in the
	 * vector
	 */
	std::vector<int> out_off_proc_rank_vec;
	/**
	 * @brief size of outgoing data from neighbors
	 */
	std::vector<int> out_off_proc_rank_size_vec;
	/**
	 * @brief map of id to local index in outgoing vector
	 */
	std::vector<int> out_off_proc_map_vec;
	/**
	 * @brief Vector of ghost patch global_indexes. Index corresponds to the patch's local index in
	 * the ghost vector.
	 */
	std::vector<int> patch_global_index_map_vec_off_proc;
	/**
	 * @brief mpi rank
	 */
	int rank;
	/**
	 * @brief The timer
	 */
	mutable std::shared_ptr<Timer> timer;

	/**
	 * @brief Give the patches local indexes.
	 *
	 * @param local_id_set true if local indexes are set by user.
	 */
	void indexDomainsLocal(bool local_id_set = false);
	/**
	 * @brief Give the patches global indexes
	 *
	 * @param global_id_set true if global indexes are set by user.
	 */
	void indexDomainsGlobal(bool global_id_set = false);

	/**
	 * @brief Give the boundary conditions local indexes.
	 */
	void indexBCLocal();
	/**
	 * @brief Give the boundary conditions global indexes
	 */
	void indexBCGlobal();

	/**
	 * @brief Give patches new global and local indexes
	 *
	 * @param local_id_set true if local indexes are set by user
	 * @param global_id_set true if global indexes are set by user
	 */
	void reIndex(bool local_id_set, bool global_id_set)
	{
		indexDomainsLocal(local_id_set);
		if (!global_id_set) {
			indexDomainsGlobal();
		}
		indexBCLocal();
		indexBCGlobal();
	}

	public:
	/**
	 * @brief Construct a new Domain object with a given PatchInfo map.
	 *
	 * @param pinfo_map map that goes from patch id to the PatchInfo pointer
	 * @param local_id_set true if local indexes are set by user
	 * @param global_id_set true if global indexes are set by user
	 */
	Domain(std::map<int, std::shared_ptr<PatchInfo<D>>> pinfo_map, std::array<int, D> ns_in,
	       int num_ghost_cells_in, bool local_id_set = false, bool global_id_set = false)
	: ns(ns_in), num_ghost_cells(num_ghost_cells_in)
	{
		num_cells_in_patch            = 1;
		num_cells_in_patch_with_ghost = 1;
		for (size_t i = 0; i < D; i++) {
			num_cells_in_patch *= ns[i];
			num_cells_in_patch_with_ghost *= (ns[i] + 2 * num_ghost_cells);
		}

		pinfo_id_map = pinfo_map;

		int num_local_domains = pinfo_id_map.size();
		MPI_Allreduce(&num_local_domains, &global_num_patches, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		reIndex(local_id_set, global_id_set);
	}
	/**
	 * @brief Get a vector of PatchInfo pointers where index in the vector corresponds to the
	 * patch's local index
	 */
	std::vector<std::shared_ptr<PatchInfo<D>>> &getPatchInfoVector()
	{
		return pinfo_vector;
	}
	/**
	 * @brief Get a vector of PatchInfo pointers where index in the vector corresponds to the
	 * patch's local index
	 */
	std::vector<std::shared_ptr<const PatchInfo<D>>> getPatchInfoVector() const
	{
		return std::vector<std::shared_ptr<const PatchInfo<D>>>(pinfo_vector.cbegin(),
		                                                        pinfo_vector.cend());
	}
	/**
	 * @brief Get map that goes form patch's id to the PatchInfo pointer
	 */
	std::map<int, std::shared_ptr<PatchInfo<D>>> &getPatchInfoMap()
	{
		return pinfo_id_map;
	}
	/**
	 * @brief Get map that goes form patch's id to the PatchInfo pointer
	 */
	std::map<int, std::shared_ptr<const PatchInfo<D>>> getPatchInfoMap() const
	{
		return std::map<int, std::shared_ptr<const PatchInfo<D>>>(pinfo_id_map.cbegin(),
		                                                          pinfo_id_map.cend());
	}
	/**
	 * @brief Get a vector of patch ids. Index in vector corresponds to the patch's local index.
	 */
	const std::vector<int> &getIdMapVec() const
	{
		return patch_id_map_vec;
	}
	/**
	 * @brief Get a vector of patch global indexes. Index in vector corresponds to the patch's
	 * global index.
	 */
	const std::vector<int> &getGlobalIndexMapVec() const
	{
		return patch_global_index_map_vec;
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
	 * @brief Set the neumann boundary conditions
	 *
	 * @param inf the function that determines boundary conditions
	 */
	void setNeumann(IsNeumannFunc<D> inf)
	{
		for (auto &p : pinfo_id_map) {
			p.second->setNeumann(inf);
		}
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
		return pinfo_vector.size();
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
		return pinfo_id_map.size() * num_cells_in_patch;
	}
	/**
	 * @brief Get get the number of local cells (including ghost cells)
	 */
	int getNumLocalCellsWithGhost() const
	{
		return pinfo_id_map.size() * num_cells_in_patch_with_ghost;
	}
	/**
	 * @brief Get get the number of local boundary condition cells
	 */
	int getNumLocalBCCells() const
	{
		return local_num_bc * num_cells_in_patch / ns[0];
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
		for (auto &p : pinfo_id_map) {
			PatchInfo<D> &d         = *p.second;
			double        patch_vol = 1;
			for (size_t i = 0; i < D; i++) {
				patch_vol *= d.spacings[i] * d.ns[i];
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

		for (auto &p : pinfo_id_map) {
			PatchInfo<D> &d = *p.second;
			for (int c = 0; c < u->getNumComponents(); c++) {
				const LocalData<D> u_data = u->getLocalData(c, d.local_index);

				double patch_sum = 0;
				nested_loop<D>(u_data.getStart(), u_data.getEnd(),
				               [&](std::array<int, D> coord) { patch_sum += u_data[coord]; });

				for (size_t i = 0; i < D; i++) {
					patch_sum *= d.spacings[i];
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
	 * @brief Set the id of the domain
	 *
	 * @param new_id the new id of the domain
	 */
	void setId(int new_id)
	{
		id = new_id;
	}
	/**
	 * @brief Get the domain's id
	 *
	 * @return int the id
	 */
	int getId()
	{
		return id;
	}
};

template <int D> void Domain<D>::indexDomainsLocal(bool local_id_set)
{
	std::vector<int>               map_vec;
	std::vector<int>               in_off_proc_map_vec;
	std::vector<int>               out_off_proc_map_vec;
	std::set<std::tuple<int, int>> in_offs;
	std::set<std::tuple<int, int>> out_offs;
	std::map<int, int>             rev_map;
	int                            curr_i = local_id_set ? pinfo_id_map.size() : 0;
	if (local_id_set) {
		if (!pinfo_id_map.empty()) {
			std::set<int> todo;
			map_vec.resize(pinfo_id_map.size());
			for (auto &p : pinfo_id_map) {
				todo.insert(p.first);
				PatchInfo<D> &d        = *p.second;
				map_vec[d.local_index] = d.id;
			}
			std::set<int> enqueued;
			while (!todo.empty()) {
				std::deque<int> queue;
				queue.push_back(*todo.begin());
				enqueued.insert(*todo.begin());
				while (!queue.empty()) {
					int i = queue.front();
					todo.erase(i);
					queue.pop_front();
					PatchInfo<D> &d = *pinfo_id_map[i];
					rev_map[i]      = d.local_index;
					auto ranks      = d.getNbrRanks();
					auto ids        = d.getNbrIds();
					for (size_t idx = 0; idx < ranks.size(); idx++) {
						int nbr_id = ids[idx];
						int rank   = ranks[idx];
						if (pinfo_id_map.count(nbr_id)) {
							if (!enqueued.count(nbr_id)) {
								enqueued.insert(nbr_id);
								queue.push_back(nbr_id);
							}
						} else {
							in_offs.insert(std::make_tuple(rank, nbr_id));
							out_offs.insert(std::make_tuple(rank, i));
						}
					}
				}
			}
		}
	} else {
		std::set<int> offs;
		if (!pinfo_id_map.empty()) {
			std::set<int> todo;
			for (auto &p : pinfo_id_map) {
				todo.insert(p.first);
			}
			std::set<int> enqueued;
			while (!todo.empty()) {
				std::deque<int> queue;
				queue.push_back(*todo.begin());
				enqueued.insert(*todo.begin());
				while (!queue.empty()) {
					int i = queue.front();
					todo.erase(i);
					queue.pop_front();
					map_vec.push_back(i);
					PatchInfo<D> &d = *pinfo_id_map[i];
					rev_map[i]      = curr_i;
					d.local_index   = curr_i;
					curr_i++;
					auto ranks = d.getNbrRanks();
					auto ids   = d.getNbrIds();
					for (size_t idx = 0; idx < ranks.size(); idx++) {
						int nbr_id = ids[idx];
						int rank   = ranks[idx];
						if (pinfo_id_map.count(nbr_id)) {
							if (!enqueued.count(nbr_id)) {
								enqueued.insert(nbr_id);
								queue.push_back(nbr_id);
							}
						} else {
							in_offs.insert(std::make_tuple(rank, nbr_id));
							out_offs.insert(std::make_tuple(rank, i));
						}
					}
				}
			}
		}
	}
	// map off proc
	for (auto pair : in_offs) {
		int rank = std::get<0>(pair);
		int id   = std::get<1>(pair);
		if (in_off_proc_rank_vec.size() == 0 || rank != in_off_proc_rank_vec.back()) {
			in_off_proc_rank_vec.push_back(rank);
			in_off_proc_rank_size_vec.push_back(1);
		} else {
			in_off_proc_rank_size_vec.back()++;
		}
		in_off_proc_map_vec.push_back(id);
		rev_map[id] = curr_i;
		curr_i++;
	}
	for (auto pair : out_offs) {
		int rank = std::get<0>(pair);
		int id   = std::get<1>(pair);
		if (out_off_proc_rank_vec.size() == 0 || rank != out_off_proc_rank_vec.back()) {
			out_off_proc_rank_vec.push_back(rank);
			out_off_proc_rank_size_vec.push_back(1);
		} else {
			out_off_proc_rank_size_vec.back()++;
		}
		out_off_proc_map_vec.push_back(id);
	}
	// set local indexes
	pinfo_vector.resize(pinfo_id_map.size());
	for (auto &p : pinfo_id_map) {
		p.second->setLocalNeighborIndexes(rev_map);
		pinfo_vector[p.second->local_index] = p.second;
	}
	// domain_rev_map          = rev_map;
	patch_global_index_map_vec          = map_vec;
	patch_id_map_vec                    = map_vec;
	this->in_off_proc_map_vec           = in_off_proc_map_vec;
	this->out_off_proc_map_vec          = out_off_proc_map_vec;
	patch_global_index_map_vec_off_proc = in_off_proc_map_vec;
}
template <int D> void Domain<D>::indexDomainsGlobal(bool global_id_set)
{
	// global indices are going to be sequentially increasing with rank
	int local_size = pinfo_id_map.size();
	int start_i;
	MPI_Scan(&local_size, &start_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	start_i -= local_size;
	patch_global_index_map_vec.resize(local_size);
	iota(patch_global_index_map_vec.begin(), patch_global_index_map_vec.end(), start_i);

	// prepare outgoing vector of data
	std::vector<int> out_data = out_off_proc_map_vec;
	for (size_t i = 0; i < out_data.size(); i++) {
		out_data[i] = patch_global_index_map_vec[pinfo_id_map[out_data[i]]->local_index];
	}
	// send outgoing messages
	using namespace std;
	vector<MPI_Request> requests;

	// recv info
	int curr_pos = 0;
	for (size_t i = 0; i < in_off_proc_rank_vec.size(); i++) {
		int         source = in_off_proc_rank_vec[i];
		int         size   = in_off_proc_rank_size_vec[i];
		MPI_Request request;
		MPI_Irecv(&patch_global_index_map_vec_off_proc[curr_pos], size, MPI_INT, source, 0,
		          MPI_COMM_WORLD, &request);
		requests.push_back(request);
		curr_pos += size;
	}
	// send info
	curr_pos = 0;
	requests.reserve(out_off_proc_rank_vec.size() + in_off_proc_rank_vec.size());
	for (size_t i = 0; i < out_off_proc_rank_vec.size(); i++) {
		int         dest = out_off_proc_rank_vec[i];
		int         size = out_off_proc_rank_size_vec[i];
		MPI_Request request;
		MPI_Isend(&out_data[curr_pos], size, MPI_INT, dest, 0, MPI_COMM_WORLD, &request);
		requests.push_back(request);
		curr_pos += size;
	}
	// wait for all
	MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

	curr_pos = 0;
	std::map<int, int> rev_map;
	for (int i : patch_global_index_map_vec) {
		rev_map[curr_pos] = i;
		curr_pos++;
	}
	for (int i : patch_global_index_map_vec_off_proc) {
		rev_map[curr_pos] = i;
		curr_pos++;
	}
	for (auto &p : pinfo_id_map) {
		p.second->setGlobalNeighborIndexes(rev_map);
	}
}
template <int D> void Domain<D>::indexBCLocal()
{
	int curr_i = 0;
	for (auto pinfo : pinfo_vector) {
		for (Side<D> s : Side<D>::getValues()) {
			if (!pinfo->hasNbr(s)) {
				pinfo->setBCLocalIndex(s, curr_i);
				patch_id_bc_map_vec.push_back(pinfo->id);
				curr_i++;
			}
		}
	}
	local_num_bc = curr_i;
}
template <int D> void Domain<D>::indexBCGlobal()
{
	int curr_i;
	MPI_Scan(&local_num_bc, &curr_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	curr_i -= local_num_bc;
	for (auto pinfo : pinfo_vector) {
		for (Side<D> s : Side<D>::getValues()) {
			if (!pinfo->hasNbr(s)) {
				pinfo->setBCGlobalIndex(s, curr_i);
				curr_i++;
			}
		}
	}
}

extern template class Domain<2>;
extern template class Domain<3>;
} // namespace ThunderEgg
#endif
