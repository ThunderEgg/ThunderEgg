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

#include "P4estDomGen.h"

#include <p4est_iterate.h>

using namespace std;
using namespace ThunderEgg;

struct my_data {
	int                rank;
	int                id;
	int                child_rank;
	std::array<int, 4> child_ids;
};
static void init_data(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
	my_data &data   = *(my_data *) quadrant->p.user_data;
	data.rank       = -1;
	data.id         = -1;
	data.child_rank = -1;
	data.child_ids.fill(-1);
}

std::shared_ptr<Domain<2>> P4estDomGen::getFinestDomain()
{
	return domain_list.front();
}
bool P4estDomGen::hasCoarserDomain()
{
	return curr_level >= 0;
}
/**
 * @brief p4est iterator to get number of levels
 */
static void get_num_levels(p4est_iter_volume_info_t *info, void *user_data)
{
	int &max_level  = *(int *) user_data;
	int  quad_level = info->quad->level;
	if (quad_level > max_level) {
		max_level = quad_level;
	}
}
/**
 * @brief p4est iterator to set gids
 */
static void set_ids(p4est_iter_volume_info_t *info, void *user_data)
{
	my_data *data = (my_data *) info->quad->p.user_data;
	data->id      = *(int *) user_data + info->quadid
	           + p4est_tree_array_index(info->p4est->trees, info->treeid)->quadrants_offset;
}
/**
 * Constructor
 */
P4estDomGen::P4estDomGen(p4est_t *p4est, const std::array<int, 2> &ns, int num_ghost_cells,
                         IsNeumannFunc<2> inf, BlockMapFunc bmf)
{
	this->num_ghost_cells = num_ghost_cells;
	this->ns              = ns;
	this->bmf             = bmf;
	this->inf             = inf;

	double x_start;
	double y_start;
	bmf(0, 0, 0, x_start, y_start);
	bmf(0, 1, 1, x_scale, y_scale);
	x_scale -= x_start;
	y_scale -= y_start;

	my_p4est = p4est_copy(p4est, false);
	p4est_reset_data(my_p4est, sizeof(my_data), init_data, nullptr);

	int max_level = 0;
	p4est_iterate_ext(my_p4est, nullptr, &max_level, get_num_levels, nullptr, nullptr, 0);
	int global_max_level = 0;
	MPI_Allreduce(&max_level, &global_max_level, 1, MPI_INT, MPI_MAX, my_p4est->mpicomm);

	int global_id_start;
	MPI_Scan(&my_p4est->local_num_quadrants, &global_id_start, 1, MPI_INT, MPI_SUM,
	         my_p4est->mpicomm);
	global_id_start -= my_p4est->local_num_quadrants;

	// set global ids
	p4est_iterate_ext(my_p4est, nullptr, &global_id_start, set_ids, nullptr, nullptr, 0);

	num_levels = global_max_level + 1;
	curr_level = global_max_level;

	// generate finest DC
	extractLevel();
}
P4estDomGen::~P4estDomGen()
{
	p4est_destroy(my_p4est);
}
/**
 * @brief iterator to set rank
 */
static void set_ranks(p4est_iter_volume_info_t *info, void *user_data)
{
	my_data *data = (my_data *) info->quad->p.user_data;
	data->rank    = info->p4est->mpirank;
}
struct create_domains_ctx {
	my_data *                                     ghost_data;
	std::array<int, 2>                            ns;
	std::map<int, std::shared_ptr<PatchInfo<2>>> *dmap;
	P4estDomGen::BlockMapFunc                     bmf;
	double                                        x_scale;
	double                                        y_scale;
	int                                           num_ghost_cells;
};
/**
 * @brief iterator to create domain objects
 */
static void create_domains(p4est_iter_volume_info_t *info, void *user_data)
{
	create_domains_ctx &ctx  = *(create_domains_ctx *) user_data;
	auto &              dmap = *ctx.dmap;

	// create domain object
	my_data *                 data = (my_data *) info->quad->p.user_data;
	shared_ptr<PatchInfo<2>> &ptr  = dmap[data->id];
	if (ptr == nullptr) {
		ptr.reset(new PatchInfo<2>);
	}
	PatchInfo<2> &pinfo = *ptr;

	pinfo.rank      = info->p4est->mpirank;
	pinfo.id        = data->id;
	pinfo.child_ids = data->child_ids;
	for (int i = 0; i < 4; i++) {
		if (pinfo.child_ids[i] != -1) {
			pinfo.child_ranks[i] = data->child_rank;
		}
	}
	pinfo.local_index
	= info->quadid + p4est_tree_array_index(info->p4est->trees, info->treeid)->quadrants_offset;

	// patch dimensions
	pinfo.ns              = ctx.ns;
	pinfo.num_ghost_cells = ctx.num_ghost_cells;

	// cell spacings
	double length     = (double) P4EST_QUADRANT_LEN(info->quad->level) / P4EST_ROOT_LEN;
	pinfo.spacings[0] = length * ctx.x_scale / ctx.ns[0];
	pinfo.spacings[1] = length * ctx.y_scale / ctx.ns[1];

	// coordinates in block
	double x = (double) info->quad->x / P4EST_ROOT_LEN;
	double y = (double) info->quad->y / P4EST_ROOT_LEN;
	ctx.bmf(info->treeid, x, y, pinfo.starts[0], pinfo.starts[1]);

	// set refinement level
	pinfo.refine_level = info->quad->level;
}
/**
 * @brief iterator to set neighbor info
 */
static void link_domains(p4est_iter_face_info_t *info, void *user_data)
{
	create_domains_ctx &ctx  = *(create_domains_ctx *) user_data;
	auto &              dmap = *ctx.dmap;
	if (info->sides.elem_count == 2) {
		p4est_iter_face_side_t side_info1 = ((p4est_iter_face_side_t *) info->sides.array)[0];
		p4est_iter_face_side_t side_info2 = ((p4est_iter_face_side_t *) info->sides.array)[1];

		auto link_side_to_side = [&](p4est_iter_face_side_t side_info1,
		                             p4est_iter_face_side_t side_info2) {
			Side<2> side(side_info1.face);
			if (side_info1.is_hanging) {
				// coarse nbr
				int nbr_id, nbr_rank;
				if (side_info2.is.full.is_ghost) {
					my_data &data = ctx.ghost_data[side_info2.is.full.quadid];
					nbr_id        = data.id;
					nbr_rank      = data.rank;
				} else {
					my_data *data = (my_data *) side_info2.is.full.quad->p.user_data;
					nbr_id        = data->id;
					nbr_rank      = info->p4est->mpirank;
				}
				for (Orthant<1> o : Orthant<1>::getValues()) {
					if (!side_info1.is.hanging.is_ghost[o.getIndex()]) {
						my_data *data
						= (my_data *) side_info1.is.hanging.quad[o.getIndex()]->p.user_data;
						PatchInfo<2> &pinfo = *dmap[data->id];
						pinfo.nbr_info[side.getIndex()].reset(new CoarseNbrInfo<2>(nbr_id, o));
						pinfo.getCoarseNbrInfo(side).rank = nbr_rank;
					}
				}
			} else if (side_info2.is_hanging) {
				// fine nbr
				if (!side_info1.is.full.is_ghost) {
					my_data *data = (my_data *) side_info1.is.full.quad->p.user_data;
					int      id   = data->id;

					std::array<int, 2> nbr_ids;
					int                ranks[2];
					for (int i = 0; i < 2; i++) {
						if (side_info2.is.hanging.is_ghost[i]) {
							my_data &data = ctx.ghost_data[side_info2.is.hanging.quadid[i]];
							nbr_ids[i]    = data.id;
							ranks[i]      = data.rank;
						} else {
							my_data *data = (my_data *) side_info2.is.hanging.quad[i]->p.user_data;
							nbr_ids[i]    = data->id;
							ranks[i]      = info->p4est->mpirank;
						}
					}
					PatchInfo<2> &pinfo = *dmap[id];
					pinfo.nbr_info[side.getIndex()].reset(new FineNbrInfo<2>(nbr_ids));
					pinfo.getFineNbrInfo(side).ranks[0] = ranks[0];
					pinfo.getFineNbrInfo(side).ranks[1] = ranks[1];
				}
			} else {
				// normal nbr
				if (!side_info1.is.full.is_ghost) {
					my_data *data = (my_data *) side_info1.is.full.quad->p.user_data;
					int      id   = data->id;
					int      nbr_id, nbr_rank;
					if (side_info2.is.full.is_ghost) {
						my_data &data = ctx.ghost_data[side_info2.is.full.quadid];
						nbr_id        = data.id;
						nbr_rank      = data.rank;
					} else {
						my_data *data = (my_data *) side_info2.is.full.quad->p.user_data;
						nbr_id        = data->id;
						nbr_rank      = info->p4est->mpirank;
					}
					PatchInfo<2> &pinfo = *dmap[id];
					pinfo.nbr_info[side.getIndex()].reset(new NormalNbrInfo<2>(nbr_id));
					pinfo.getNormalNbrInfo(side).rank = nbr_rank;
				}
			}
		};

		// 1 to 2
		link_side_to_side(side_info1, side_info2);
		// 2 to 1
		link_side_to_side(side_info2, side_info1);
	}
}

struct coarsen_domains_ctx {
	int                                 level;
	map<int, shared_ptr<PatchInfo<2>>> *dmap;
};
/**
 * @brief iterator to coarsen to next level
 */
static int coarsen(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quads[])
{
	coarsen_domains_ctx &ctx = *(coarsen_domains_ctx *) p4est->user_pointer;
	if (quads[1] == nullptr) {
		my_data *data      = (my_data *) quads[0]->p.user_data;
		data->child_ids[0] = data->id;
		data->child_rank   = data->rank;
		return false;
	} else if (quads[0]->level > ctx.level) {
		return true;
	} else {
		for (int i = 0; i < 4; i++) {
			my_data *data      = (my_data *) quads[i]->p.user_data;
			data->child_ids[0] = data->id;
			data->child_rank   = data->rank;
		}
		return false;
	}
}

/**
 * @brief set parent info in finer level
 */
static void coarsen_replace(p4est_t *p4est, p4est_topidx_t which_tree, int num_outgoing,
                            p4est_quadrant_t *outgoing[], int num_incoming,
                            p4est_quadrant_t *incoming[])
{
	coarsen_domains_ctx &ctx       = *(coarsen_domains_ctx *) p4est->user_pointer;
	auto &               dmap      = *ctx.dmap;
	my_data *            data      = (my_data *) incoming[0]->p.user_data;
	my_data *            fine_data = (my_data *) outgoing[0]->p.user_data;
	data->id                       = fine_data->id;
	data->child_rank               = fine_data->rank;

	for (Orthant<2> o : Orthant<2>::getValues()) {
		my_data *     fine_data       = (my_data *) outgoing[o.getIndex()]->p.user_data;
		PatchInfo<2> &pinfo           = *dmap.at(fine_data->id);
		pinfo.orth_on_parent          = o;
		pinfo.parent_id               = data->id;
		data->child_ids[o.getIndex()] = fine_data->id;
	}
}
void P4estDomGen::extractLevel()
{
	// coarsen previous
	if (domain_list.size() > 0) {
		coarsen_domains_ctx ctx = {curr_level, &domain_list.back()->getPatchInfoMap()};
		my_p4est->user_pointer  = &ctx;
		p4est_coarsen_ext(my_p4est, false, true, coarsen, init_data, coarsen_replace);
		for (auto p : domain_list.back()->getPatchInfoMap()) {
			PatchInfo<2> &pinfo = *p.second;
			if (!pinfo.hasCoarseParent()) {
				pinfo.parent_id = pinfo.id;
			}
		}
		p4est_partition_ext(my_p4est, true, nullptr);
	}
	// set ranks
	p4est_iterate_ext(my_p4est, nullptr, nullptr, set_ranks, nullptr, nullptr, 0);

	// get ids of ghosts
	p4est_ghost_t *ghost = p4est_ghost_new(my_p4est, P4EST_CONNECT_FACE);
	my_data        ghost_data[ghost->ghosts.elem_count];
	p4est_ghost_exchange_data(my_p4est, ghost, ghost_data);

	std::map<int, std::shared_ptr<PatchInfo<2>>> new_level;
	create_domains_ctx ctx = {ghost_data, ns, &new_level, bmf, x_scale, y_scale, num_ghost_cells};
	// create domain objects and set neighbor information
	p4est_iterate_ext(my_p4est, ghost, &ctx, create_domains, link_domains, nullptr, 0);

	p4est_ghost_destroy(ghost);

	for (auto p : new_level) {
		p.second->setPtrs(new_level);
	}
	if (domain_list.size() > 0) {
		auto &old_level = domain_list.back()->getPatchInfoMap();
		// update parent ranks
		std::map<int, int> id_rank_map;
		// get outgoing information
		std::map<int, std::set<int>> out_info;
		for (auto pair : new_level) {
			id_rank_map[pair.first] = my_p4est->mpirank;
			auto &pinfo             = pair.second;
			for (int i = 0; i < 4; i++) {
				if (pinfo->child_ranks[i] != -1 && pinfo->child_ranks[i] != my_p4est->mpirank) {
					out_info[pinfo->child_ranks[i]].insert(pinfo->id);
				}
			}
		}
		// get incoming information
		std::set<int> incoming_ids;
		for (auto &pair : old_level) {
			if (!new_level.count(pair.second->parent_id)) {
				incoming_ids.insert(pair.second->parent_id);
			}
		}
		// send info
		std::deque<std::vector<int>> buffers;
		std::vector<MPI_Request>     send_requests;
		for (auto &pair : out_info) {
			int              dest = pair.first;
			std::vector<int> buffer(pair.second.begin(), pair.second.end());
			buffers.push_back(buffer);
			MPI_Request request;
			MPI_Isend(buffer.data(), buffer.size(), MPI_INT, dest, 0, MPI_COMM_WORLD, &request);
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
		MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
		// update rank info
		for (auto &pair : old_level) {
			pair.second->parent_rank = id_rank_map.at(pair.second->parent_id);
		}
	}
	// create Domain object
	domain_list.push_back(
	shared_ptr<Domain<2>>(new Domain<2>(new_level, ns, num_ghost_cells, true)));
	domain_list.back()->setNeumann(inf);

	curr_level--;
}
std::shared_ptr<Domain<2>> P4estDomGen::getCoarserDomain()
{
	if (curr_level >= 0) {
		extractLevel();
		return domain_list.back();
	} else {
		return nullptr;
	}
}
