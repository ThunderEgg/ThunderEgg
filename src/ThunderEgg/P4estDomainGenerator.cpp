/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

#include "P4estDomainGenerator.h"

#include <p4est_iterate.h>

using namespace std;
using namespace ThunderEgg;

namespace
{
struct IterateFunctions {
	std::function<void(p4est_iter_volume_info_t *)> iter_volume;
	std::function<void(p4est_iter_face_info_t *)>   iter_face;
	std::function<void(p4est_iter_corner_info_t *)> iter_corner;
};

void IterVolumeWrap(p4est_iter_volume_info_t *info, void *user_data)
{
	const IterateFunctions *functions = (IterateFunctions *) user_data;
	functions->iter_volume(info);
}

void IterFaceWrap(p4est_iter_face_info_t *info, void *user_data)
{
	const IterateFunctions *functions = (IterateFunctions *) user_data;
	functions->iter_face(info);
}

void IterCornerWrap(p4est_iter_corner_info_t *info, void *user_data)
{
	const IterateFunctions *functions = (IterateFunctions *) user_data;
	functions->iter_corner(info);
}

void p4est_iterate_ext(p4est_t *                                                    p4est,
                       p4est_ghost_t *                                              ghost_layer,
                       const std::function<void(const p4est_iter_volume_info_t *)> &iter_volume,
                       const std::function<void(const p4est_iter_face_info_t *)> &  iter_face,
                       const std::function<void(const p4est_iter_corner_info_t *)> &iter_corner,
                       int                                                          remote)
{
	IterateFunctions functions = {iter_volume, iter_face, iter_corner};
	p4est_iterate_ext(p4est,
	                  ghost_layer,
	                  &functions,
	                  iter_volume ? &IterVolumeWrap : nullptr,
	                  iter_face ? &IterFaceWrap : nullptr,
	                  iter_corner ? &IterCornerWrap : nullptr,
	                  remote);
}

struct CoarsenFunctions {
	std::function<int(p4est_topidx_t, p4est_quadrant_t *[])>                                  coarsen;
	std::function<void(p4est_topidx_t, int, p4est_quadrant_t *[], int, p4est_quadrant_t *[])> replace;
};

int CoarsenWrap(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrants[])
{
	const CoarsenFunctions *functions = (CoarsenFunctions *) p4est->user_pointer;
	return functions->coarsen(which_tree, quadrants);
}

void ReplaceWrap(p4est_t *         p4est,
                 p4est_topidx_t    which_tree,
                 int               num_outgoing,
                 p4est_quadrant_t *outgoing[],
                 int               num_incoming,
                 p4est_quadrant_t *incoming[])
{
	const CoarsenFunctions *functions = (CoarsenFunctions *) p4est->user_pointer;
	functions->replace(which_tree, num_outgoing, outgoing, num_incoming, incoming);
}

void p4est_coarsen_ext(p4est_t *                                                                                        p4est,
                       int                                                                                              coarsen_recursive,
                       int                                                                                              callback_orphans,
                       const std::function<int(p4est_topidx_t, p4est_quadrant_t *[])> &                                 coarsen,
                       p4est_init_t                                                                                     init_fn,
                       const std::function<void(p4est_topidx_t, int, p4est_quadrant_t *[], int, p4est_quadrant_t *[])> &replace)
{
	CoarsenFunctions functions = {coarsen, replace};
	p4est->user_pointer        = &functions;
	p4est_coarsen_ext(p4est, coarsen_recursive, callback_orphans, coarsen ? &CoarsenWrap : nullptr, init_fn, replace ? &ReplaceWrap : nullptr);
	p4est->user_pointer = nullptr;
}

int GetMaxLevel(p4est_t *p4est)
{
	int max_level = 0;

	p4est_iterate_ext(
	p4est, nullptr, [&](const p4est_iter_volume_info_t *info) { max_level = max(max_level, (int) info->quad->level); }, nullptr, nullptr, false);

	int global_max_level;

	MPI_Allreduce(&max_level, &global_max_level, 1, MPI_INT, MPI_MAX, p4est->mpicomm);
	return global_max_level;
}

struct Data {
	int                id;
	int                rank;
	std::array<int, 4> child_ids;
	std::array<int, 4> child_ranks;
	PatchInfo<2> *     pinfo;
};

void InitData(p4est_t *, p4est_topidx_t, p4est_quadrant_t *quadrant)
{
	Data &data = *(Data *) quadrant->p.user_data;
	data.id    = -1;
	data.rank  = -1;
	data.child_ids.fill(-1);
	data.child_ranks.fill(-1);
	data.pinfo = nullptr;
}

void SetRanks(p4est_t *p4est)
{
	int rank;
	MPI_Comm_rank(p4est->mpicomm, &rank);

	p4est_iterate_ext(
	p4est,
	nullptr,
	[=](const p4est_iter_volume_info_t *info) {
		Data *data = (Data *) info->quad->p.user_data;
		data->rank = rank;
	},
	nullptr,
	nullptr,
	false);
}

void SetIds(p4est_t *p4est)
{
	int curr_global_id;
	MPI_Scan(&p4est->local_num_quadrants, &curr_global_id, 1, MPI_INT, MPI_SUM, p4est->mpicomm);
	curr_global_id -= p4est->local_num_quadrants;

	p4est_iterate_ext(
	p4est,
	nullptr,
	[&](const p4est_iter_volume_info_t *info) {
		Data *data = (Data *) info->quad->p.user_data;
		data->id   = curr_global_id;
		curr_global_id++;
	},
	nullptr,
	nullptr,
	false);
}

} // namespace

/*
 * constructor
 */
P4estDomainGenerator::P4estDomainGenerator(p4est_t *p4est, const std::array<int, 2> &ns, int num_ghost_cells, const BlockMapFunc &bmf)
: ns(ns),
  num_ghost_cells(num_ghost_cells),
  bmf(bmf)
{
	my_p4est = p4est_copy(p4est, false);

	p4est_reset_data(my_p4est, sizeof(Data), InitData, nullptr);

	curr_level = GetMaxLevel(my_p4est);

	SetIds(my_p4est);

	extractLevel();
}

P4estDomainGenerator::~P4estDomainGenerator()
{
	p4est_destroy(my_p4est);
}

P4estDomainGenerator::P4estDomainGenerator(const P4estDomainGenerator &other)
: domain_patches(other.domain_patches),
  ns(other.ns),
  num_ghost_cells(other.num_ghost_cells),
  curr_level(other.curr_level),
  bmf(other.bmf)
{
	my_p4est = p4est_copy(other.my_p4est, true);
}

P4estDomainGenerator &P4estDomainGenerator::operator=(const P4estDomainGenerator &other)
{
	domain_patches  = other.domain_patches;
	ns              = other.ns;
	num_ghost_cells = other.num_ghost_cells;
	curr_level      = other.curr_level;
	bmf             = other.bmf;
	p4est_destroy(my_p4est);
	my_p4est = p4est_copy(other.my_p4est, true);
	return *this;
}
void P4estDomainGenerator::extractLevel()
{
	if (domain_patches.size() > 0) {
		coarsenTree();
	}

	SetRanks(my_p4est);

	createPatchInfos();

	linkNeighbors();

	if (domain_patches.size() == 2) {
		updateParentRanksOfPreviousDomain();
	}

	curr_level--;
}
void P4estDomainGenerator::coarsenTree()
{
	p4est_coarsen_ext(
	my_p4est,
	false,
	true,
	[&](p4est_topidx_t, p4est_quadrant_t *quadrant[]) {
		if (quadrant[1] == nullptr) {
			// update child data to point to self
			Data *data                  = (Data *) quadrant[0]->p.user_data;
			data->child_ids[0]          = data->id;
			data->child_ranks[0]        = data->rank;
			data->pinfo->parent_id      = data->id;
			data->pinfo->orth_on_parent = Orthant<2>::null();
			return false;
		} else if (quadrant[0]->level > curr_level) {
			// update child datas to point to parent
			const Data *bsw_data = (Data *) quadrant[0]->p.user_data;
			for (unsigned char i = 0; i < 4; i++) {
				Data *data                  = (Data *) quadrant[i]->p.user_data;
				data->child_ids[0]          = data->id;
				data->child_ranks[0]        = data->rank;
				data->pinfo->parent_id      = bsw_data->id;
				data->pinfo->orth_on_parent = Orthant<2>(i);
			}
			return true;
		} else {
			// update child datas to point to self
			for (int i = 0; i < 4; i++) {
				Data *data                  = (Data *) quadrant[i]->p.user_data;
				data->child_ids[0]          = data->id;
				data->child_ranks[0]        = data->rank;
				data->pinfo->parent_id      = data->id;
				data->pinfo->orth_on_parent = Orthant<2>::null();
			}
			return false;
		}
	},
	InitData,
	[](p4est_topidx_t, int, p4est_quadrant_t *outgoing[], int, p4est_quadrant_t *incoming[]) {
		Data *      data         = (Data *) incoming[0]->p.user_data;
		const Data *fine_sw_data = (Data *) outgoing[0]->p.user_data;
		data->id                 = fine_sw_data->id;
		for (int i = 0; i < 4; i++) {
			const Data *fine_data = (Data *) outgoing[i]->p.user_data;
			data->child_ids[i]    = fine_data->id;
			data->child_ranks[i]  = fine_data->rank;
		}
	});
	p4est_partition_ext(my_p4est, true, nullptr);
}

void P4estDomainGenerator::createPatchInfos()
{
	domain_patches.emplace_front(my_p4est->local_num_quadrants);
	p4est_iterate_ext(
	my_p4est,
	nullptr,
	[&](const p4est_iter_volume_info_t *info) {
		Data *data  = (Data *) info->quad->p.user_data;
		data->pinfo = &domain_patches.front()[info->quadid + p4est_tree_array_index(info->p4est->trees, info->treeid)->quadrants_offset];

		data->pinfo->rank            = info->p4est->mpirank;
		data->pinfo->id              = data->id;
		data->pinfo->ns              = ns;
		data->pinfo->num_ghost_cells = num_ghost_cells;
		data->pinfo->refine_level    = info->quad->level;

		data->pinfo->child_ids   = data->child_ids;
		data->pinfo->child_ranks = data->child_ranks;

		double unit_x = (double) info->quad->x / P4EST_ROOT_LEN;
		double unit_y = (double) info->quad->y / P4EST_ROOT_LEN;

		bmf(info->treeid, unit_x, unit_y, data->pinfo->starts[0], data->pinfo->starts[1]);

		double upper_unit_x = (double) (info->quad->x + P4EST_QUADRANT_LEN(info->quad->level)) / P4EST_ROOT_LEN;
		double upper_unit_y = (double) (info->quad->y + P4EST_QUADRANT_LEN(info->quad->level)) / P4EST_ROOT_LEN;

		bmf(info->treeid, upper_unit_x, upper_unit_y, data->pinfo->spacings[0], data->pinfo->spacings[1]);

		for (int i = 0; i < 2; i++) {
			data->pinfo->spacings[i] -= data->pinfo->starts[i];
			data->pinfo->spacings[i] /= ns[i];
		}
	},
	nullptr,
	nullptr,
	false);
}

/*
 * Functions for linkNeighbors
 */
namespace
{
/**
 * @brief set the CoarseNbrInfos on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data
 */
void addCoarseNbrInfo(const p4est_iter_face_side_t &side1, const p4est_iter_face_side_t &side2, const Data *ghost_data)
{
	const Data *nbr_data;
	if (side2.is.full.is_ghost) {
		nbr_data = ghost_data + side2.is.full.quadid;
	} else {
		nbr_data = (Data *) side2.is.full.quad->p.user_data;
	}
	for (unsigned char i = 0; i < 2; i++) {
		if (!side1.is.hanging.is_ghost[i]) {
			Data *data = (Data *) side1.is.hanging.quad[i]->p.user_data;

			Side<2> side(side1.face);

			CoarseNbrInfo<1> *nbr_info = new CoarseNbrInfo<1>();
			nbr_info->id               = nbr_data->id;
			nbr_info->rank             = nbr_data->rank;
			nbr_info->orth_on_coarse   = Orthant<1>(i);
			data->pinfo->setNbrInfo(side, nbr_info);
		}
	}
}

/**
 * @brief set the FineNbrInfo on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data
 */
void addFineNbrInfo(const p4est_iter_face_side_t &side1, const p4est_iter_face_side_t &side2, const Data *ghost_data)
{
	if (!side1.is.full.is_ghost) {
		Data *      data = (Data *) side1.is.full.quad->p.user_data;
		const Data *nbr_datas[2];
		for (int i = 0; i < 2; i++) {
			if (side2.is.hanging.is_ghost[i]) {
				nbr_datas[i] = ghost_data + side2.is.hanging.quadid[i];
			} else {
				nbr_datas[i] = (Data *) side2.is.hanging.quad[i]->p.user_data;
			}
		}

		Side<2> side(side1.face);

		FineNbrInfo<1> *nbr_info = new FineNbrInfo<1>();
		for (int i = 0; i < 2; i++) {
			nbr_info->ids[i]   = nbr_datas[i]->id;
			nbr_info->ranks[i] = nbr_datas[i]->rank;
		}
		data->pinfo->setNbrInfo(side, nbr_info);
	}
}

/**
 * @brief Set the NormalNbrInfo on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data
 */
void addNormalNbrInfo(const p4est_iter_face_side_t &side1, const p4est_iter_face_side_t &side2, const Data *ghost_data)
{
	if (!side1.is.full.is_ghost) {
		Data *      data = (Data *) side1.is.full.quad->p.user_data;
		const Data *nbr_data;
		if (side2.is.full.is_ghost) {
			nbr_data = ghost_data + side2.is.full.quadid;
		} else {
			nbr_data = (Data *) side2.is.full.quad->p.user_data;
		}

		Side<2> side(side1.face);

		NormalNbrInfo<1> *nbr_info = new NormalNbrInfo<1>(nbr_data->id);
		nbr_info->rank             = nbr_data->rank;
		data->pinfo->setNbrInfo(side, nbr_info);
	}
}

/**
 * @brief Add information from side2 to NbrInfo objects on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data p4est ghost data array
 */
void addNbrInfo(const p4est_iter_face_side_t &side1, const p4est_iter_face_side_t &side2, const Data *ghost_data)
{
	if (side1.is_hanging) {
		addCoarseNbrInfo(side1, side2, ghost_data);
	} else if (side2.is_hanging) {
		addFineNbrInfo(side1, side2, ghost_data);
	} else {
		addNormalNbrInfo(side1, side2, ghost_data);
	}
}

/**
 * @brief set the CornerFineNbrInfo on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data
 */
void addCornerFineNbrInfo(const p4est_iter_corner_side_t &side1, const p4est_iter_corner_side_t &side2, const Data *ghost_data)
{
	if (!side1.is_ghost) {
		Data *      data = (Data *) side1.quad->p.user_data;
		const Data *nbr_data;
		if (side2.is_ghost) {
			nbr_data = ghost_data + side2.quadid;
		} else {
			nbr_data = (Data *) side2.quad->p.user_data;
		}

		Corner<2> corner(side1.corner);

		FineNbrInfo<0> *nbr_info = new FineNbrInfo<0>();
		nbr_info->ids[0]         = nbr_data->id;
		nbr_info->ranks[0]       = nbr_data->rank;
		data->pinfo->setNbrInfo(corner, nbr_info);
	}
}

/**
 * @brief set the CornerCoarseNbrInfos on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data
 */
void addCornerCoarseNbrInfo(const p4est_iter_corner_side_t &side1, const p4est_iter_corner_side_t &side2, const Data *ghost_data)
{
	if (!side1.is_ghost) {
		const Data *nbr_data;
		if (side2.is_ghost) {
			nbr_data = ghost_data + side2.quadid;
		} else {
			nbr_data = (Data *) side2.quad->p.user_data;
		}
		Data *data = (Data *) side1.quad->p.user_data;

		Corner<2> corner(side1.corner);

		CoarseNbrInfo<0> *nbr_info = new CoarseNbrInfo<0>();
		nbr_info->id               = nbr_data->id;
		nbr_info->rank             = nbr_data->rank;
		nbr_info->orth_on_coarse   = Orthant<0>(0);
		data->pinfo->setNbrInfo(corner, nbr_info);
	}
}

/**
 * @brief Set the CornerNormalNbrInfo on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data
 */
void addCornerNormalNbrInfo(const p4est_iter_corner_side_t &side1, const p4est_iter_corner_side_t &side2, const Data *ghost_data)
{
	if (!side1.is_ghost) {
		Data *      data = (Data *) side1.quad->p.user_data;
		const Data *nbr_data;
		if (side2.is_ghost) {
			nbr_data = ghost_data + side2.quadid;
		} else {
			nbr_data = (Data *) side2.quad->p.user_data;
		}

		Corner<2> corner(side1.corner);

		NormalNbrInfo<0> *nbr_info = new NormalNbrInfo<0>(nbr_data->id);
		nbr_info->rank             = nbr_data->rank;
		data->pinfo->setNbrInfo(corner, nbr_info);
	}
}

/**
 * @brief Add information from side2 to NbrInfo objects on side1
 *
 * @param side1
 * @param side2
 * @param ghost_data p4est ghost data array
 */
void addCornerNbrInfo(const p4est_iter_corner_side_t &side1, const p4est_iter_corner_side_t &side2, const Data *ghost_data)
{
	if (side1.quad->level > side2.quad->level) {
		addCornerCoarseNbrInfo(side1, side2, ghost_data);
	} else if (side1.quad->level < side2.quad->level) {
		addCornerFineNbrInfo(side1, side2, ghost_data);
	} else {
		addCornerNormalNbrInfo(side1, side2, ghost_data);
	}
}
} // namespace

void P4estDomainGenerator::linkNeighbors()
{
	p4est_ghost_t *ghost = p4est_ghost_new(my_p4est, P4EST_CONNECT_CORNER);

	Data ghost_data[ghost->ghosts.elem_count];
	p4est_ghost_exchange_data(my_p4est, ghost, ghost_data);

	p4est_iterate_ext(
	my_p4est,
	ghost,
	nullptr,
	[&](const p4est_iter_face_info_t *info) {
		if (info->sides.elem_count == 2) {
			p4est_iter_face_side_t side1 = ((p4est_iter_face_side_t *) info->sides.array)[0];
			p4est_iter_face_side_t side2 = ((p4est_iter_face_side_t *) info->sides.array)[1];
			addNbrInfo(side1, side2, ghost_data);
			addNbrInfo(side2, side1, ghost_data);
		}
	},
	[&](const p4est_iter_corner_info *info) {
		if (info->sides.elem_count == 4) {
			for (Corner<2> c : Corner<2>::getValues()) {
				p4est_iter_corner_side_t side1 = ((p4est_iter_corner_side_t *) info->sides.array)[c.getIndex()];
				p4est_iter_corner_side_t side2 = ((p4est_iter_corner_side_t *) info->sides.array)[c.opposite().getIndex()];
				addCornerNbrInfo(side1, side2, ghost_data);
			}
		}
	},
	false);

	p4est_ghost_destroy(ghost);
}

void P4estDomainGenerator::updateParentRanksOfPreviousDomain()
{
	std::vector<PatchInfo<2>> &      old_level = domain_patches.back();
	const std::vector<PatchInfo<2>> &new_level = domain_patches.front();

	std::set<int> local_new_level_ids;
	for (const PatchInfo<2> &pinfo : new_level) {
		local_new_level_ids.insert(pinfo.id);
	}

	// update parent ranks
	std::map<int, int> id_rank_map;
	// get outgoing information
	std::map<int, std::set<int>> out_info;
	for (const PatchInfo<2> &pinfo : new_level) {
		id_rank_map[pinfo.id] = my_p4est->mpirank;
		for (int i = 0; i < 4; i++) {
			if (pinfo.child_ranks[i] != -1 && pinfo.child_ranks[i] != my_p4est->mpirank) {
				out_info[pinfo.child_ranks[i]].insert(pinfo.id);
			}
		}
	}
	// get incoming information
	std::set<int> incoming_ids;
	for (const PatchInfo<2> &pinfo : old_level) {
		if (!local_new_level_ids.count(pinfo.parent_id)) {
			incoming_ids.insert(pinfo.parent_id);
		}
	}
	// send info
	std::deque<std::vector<int>> buffers;
	std::vector<MPI_Request>     send_requests;
	for (const auto &pair : out_info) {
		int              dest = pair.first;
		std::vector<int> buffer(pair.second.begin(), pair.second.end());
		buffers.push_back(buffer);
		MPI_Request request;
		MPI_Isend(buffer.data(), (int) buffer.size(), MPI_INT, dest, 0, MPI_COMM_WORLD, &request);
		send_requests.push_back(request);
	}
	// recv info
	while (incoming_ids.size()) {
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		int size;
		MPI_Get_count(&status, MPI_INT, &size);
		int *buffer = new int[size];

		MPI_Recv(buffer, size, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);

		for (int i = 0; i < size; i++) {
			id_rank_map[buffer[i]] = status.MPI_SOURCE;
			incoming_ids.erase(buffer[i]);
		}

		delete[] buffer;
	}
	// wait for all
	MPI_Waitall((int) send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
	// update rank info
	for (PatchInfo<2> &pinfo : old_level) {
		pinfo.parent_rank = id_rank_map.at(pinfo.parent_id);
	}
}

Domain<2> P4estDomainGenerator::getCoarserDomain()
{
	if (curr_level >= 0) {
		extractLevel();
	}
	Communicator comm(my_p4est->mpicomm);
	Domain<2>    domain(comm, id, ns, num_ghost_cells, domain_patches.back().begin(), domain_patches.back().end());
	domain_patches.pop_back();
	id++;
	return domain;
}

Domain<2> P4estDomainGenerator::getFinestDomain()
{
	if (curr_level >= 0) {
		extractLevel();
	}
	Communicator comm(my_p4est->mpicomm);
	Domain<2>    domain(comm, id, ns, num_ghost_cells, domain_patches.back().begin(), domain_patches.back().end());
	domain_patches.pop_back();
	id++;
	return domain;
}

bool P4estDomainGenerator::hasCoarserDomain()
{
	return curr_level > 0;
}
