/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_GMG_INTERLEVELCOMM_H
#define THUNDEREGG_GMG_INTERLEVELCOMM_H

#include <Thunderegg/Domain.h>
#include <Thunderegg/ValVector.h>

namespace Thunderegg
{
namespace GMG
{
/**
 * @brief Exception that the InterLevelComm class trows
 */
struct InterLevelCommException : std::runtime_error {
	InterLevelCommException(std::string message) : std::runtime_error(message){};
};
/**
 * @brief Creates a mapping from fine to coarse levels.
 */
template <size_t D> class InterLevelComm
{
	private:
	std::array<int, D>                                               ns;
	int                                                              num_ghost_cells;
	std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> patches_with_local_parent;
	std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> patches_with_ghost_parent;
	std::vector<int>                                                 ghost_vector_rank_vector;
	std::vector<std::pair<int, std::vector<int>>>                    vector_rank_local_indexes;
	std::vector<std::pair<int, std::vector<int>>>                    ghost_rank_local_indexes;
	bool                                                             communicating = false;
	bool                                                             sending       = false;

	std::shared_ptr<const Vector<D>> current_vector;
	std::shared_ptr<const Vector<D>> current_ghost_vector;

	public:
	/**
	 * @brief Create a new InterLevelComm object.
	 *
	 * @param coarse_domain the coarser DomainCollection.
	 * @param fine_domain the finer DomainCollection.
	 */
	InterLevelComm(std::shared_ptr<const Domain<D>> coarser_domain,
	               std::shared_ptr<const Domain<D>> finer_domain)
	: ns(finer_domain->getNs()), num_ghost_cells(finer_domain->getNumGhostCells())
	{
		// sort into patches with local parents and patches with ghost parents
		std::deque<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> local_parents;
		std::deque<std::shared_ptr<const PatchInfo<D>>>                 ghost_parents;
		std::set<int>                                                   ghost_parents_ids;

		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		// TODO this has to be changed when domain class is updated
		auto cd_map = coarser_domain->getPatchInfoMap();
		for (auto patch : finer_domain->getPatchInfoVector()) {
			if (patch->parent_rank == rank) {
				local_parents.emplace_back(cd_map.at(patch->parent_id)->local_index, patch);
			} else {
				ghost_parents.push_back(patch);
				ghost_parents_ids.insert(patch->parent_id);
			}
		}

		// fill in local vector
		patches_with_local_parent
		= std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>>(local_parents.begin(),
		                                                                   local_parents.end());

		// fill in ghost vector
		// first, assign local indexes in ghost vector for ghost patches
		std::map<int, int> id_ghost_vector_local_index_map;
		int                index = 0;
		for (int id : ghost_parents_ids) {
			id_ghost_vector_local_index_map[id] = index;
			index++;
		}

		patches_with_ghost_parent.reserve(ghost_parents.size());
		ghost_vector_rank_vector.resize(ghost_parents_ids.size());
		for (auto patch : ghost_parents) {
			int ghost_local_index = id_ghost_vector_local_index_map[patch->parent_id];

			patches_with_ghost_parent.emplace_back(ghost_local_index, patch);

			ghost_vector_rank_vector[ghost_local_index] = patch->parent_rank;
		}

		// find local coarse patches that are ghost paches on other ranks
		std::map<int, std::set<std::pair<int, int>>> ranks_and_local_patches;
		for (auto pinfo : coarser_domain->getPatchInfoVector()) {
			for (int child_rank : pinfo->child_ranks) {
				if (child_rank != -1 && child_rank != rank)
					ranks_and_local_patches[child_rank].emplace(pinfo->id, pinfo->local_index);
			}
		}

		// the map sould have sorted patches by id, the other processors will expect things in this
		// order
		vector_rank_local_indexes.reserve(ranks_and_local_patches.size());
		for (auto pair : ranks_and_local_patches) {
			std::vector<int> local_indexes;
			local_indexes.reserve(pair.second.size());
			for (auto id_local_index_pair : pair.second) {
				local_indexes.push_back(id_local_index_pair.second);
			}
			vector_rank_local_indexes.emplace_back(pair.first, local_indexes);
		}
	}

	/**
	 * @brief Allocate a new vector for ghost patch values
	 *
	 * @return the newly allocated vector.
	 */
	std::shared_ptr<Vector<D>> getNewGhostVector() const
	{
		return std::make_shared<ValVector<D>>(ns, num_ghost_cells, ghost_vector_rank_vector.size());
	}

	const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &
	getPatchesWithLocalParent() const
	{
		return patches_with_local_parent;
	}
	const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &
	getPatchesWithGhostParent() const
	{
		return patches_with_ghost_parent;
	}

	/**
	 * @brief Get a PETSc VecScatter object that will scatter values from the coarse vector to the
	 * processes that need them for interpolation / restriction.
	 *
	 * @return the VecScatter object
	 */
	void sendGhostPatchesStart(std::shared_ptr<Vector<D>>       vector,
	                           std::shared_ptr<const Vector<D>> ghost_vector)
	{
		if (communicating) {
			if (sending) {
				throw InterLevelCommException(
				"InterLevelComm has a sendGhostPatches posted that is unfinished");
			} else {
				throw InterLevelCommException(
				"InterLevelComm has a getGhostPatches posted that is unfinished");
			}
		}

		// keep track of what vectors we are using
		current_ghost_vector = ghost_vector;
		current_vector       = vector;

		// set state
		communicating = true;
		sending       = true;
	}
	void sendGhostPatchesFinish(std::shared_ptr<Vector<D>>       vector,
	                            std::shared_ptr<const Vector<D>> ghost_vector)
	{
		if (!communicating) {
			throw InterLevelCommException(
			"InterLevelComm cannot finish sendGhostPatches since communication was not started");
		} else if (!sending) {
			throw InterLevelCommException(
			"InterLevelComm sendGhostPatchesFinish is being called after getGhostPatchesStart was called");
		}
		if (vector != current_vector) {
			throw InterLevelCommException(
			"InterLevelComm sendGhostPatchesFinish is being called with a different vector than when sendGhostPatchesStart was called");
		}
		if (ghost_vector != current_ghost_vector) {
			throw InterLevelCommException(
			"InterLevelComm senGhostPatchesFinish is being called with a different ghost vector than when sendGhostPatchesStart was called");
		}

		// set state
		communicating        = false;
		current_ghost_vector = nullptr;
		current_vector       = nullptr;
	}
	void getGhostPatchesStart(std::shared_ptr<const Vector<D>> vector,
	                          std::shared_ptr<Vector<D>>       ghost_vector)
	{
		if (communicating) {
			if (sending) {
				throw InterLevelCommException(
				"InterLevelComm has a sendGhostPatches posted that is unfinished");
			} else {
				throw InterLevelCommException(
				"InterLevelComm has a getGhostPatches posted that is unfinished");
			}
		}

		// keep track of what vectors we are using
		current_ghost_vector = ghost_vector;
		current_vector       = vector;

		// set state
		communicating = true;
		sending       = false;
	}
	void getGhostPatchesFinish(std::shared_ptr<const Vector<D>> vector,
	                           std::shared_ptr<Vector<D>>       ghost_vector)
	{
		if (!communicating) {
			throw InterLevelCommException(
			"InterLevelComm cannot finish sendGhostPatches since communication was not started");
		} else if (sending) {
			throw InterLevelCommException(
			"InterLevelComm getGhostPatchesFinish is being called after sendGhostPatchesStart was called");
		}
		if (vector != current_vector) {
			throw InterLevelCommException(
			"InterLevelComm getGhostPatchesFinish is being called with a different vector than when getGhostPatchesStart was called");
		}
		if (ghost_vector != current_ghost_vector) {
			throw InterLevelCommException(
			"InterLevelComm getGhostPatchesFinish is being called with a different ghost vector than when getGhostPatchesStart was called");
		}

		// set state
		communicating        = false;
		current_ghost_vector = nullptr;
		current_vector       = nullptr;
	}
}; // namespace GMG
} // namespace GMG
} // namespace Thunderegg
#endif