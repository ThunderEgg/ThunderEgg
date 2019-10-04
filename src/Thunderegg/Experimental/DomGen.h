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

#ifndef THUNDEREGG_EXPERIMENTAL_DOMGEN_H
#define THUNDEREGG_EXPERIMENTAL_DOMGEN_H
#include <Thunderegg/DomainGenerator.h>
#include <Thunderegg/Experimental/OctTree.h>
#include <list>
#include <zoltan.h>
namespace Thunderegg
{
namespace Experimental
{
/**
 * @brief Creates domains from a Thunderegg Tree
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class DomGen : public DomainGenerator<D>
{
	private:
	/**
	 * @brief number of ghost cells on each side of the patch
	 */
	int num_ghost_cells;
	/**
	 * @brief Whether Neumann boundary conditions are being used
	 */
	bool neumann;
	/**
	 * @brief finest is stored in front
	 */
	std::list<std::shared_ptr<Domain<D>>> domain_list;
	/**
	 * @brief the level of the generated domain. starts with num_levels-1
	 */
	int curr_level;
	/**
	 * @brief number of levels in the tree
	 */
	int num_levels;
	/**
	 * @brief the tree
	 */
	Tree<D> t;
	/**
	 * @brief map of id to rank
	 */
	std::map<int, int> id_rank_map;
	/**
	 * @brief the number of cells in each direction.
	 */
	std::array<int, D> ns;
	/**
	 * @brief extract a new coarser level
	 */
	void extractLevel();
	/**
	 * @brief maps patch id to PatchInfo pointer
	 */
	using PInfoMap = std::map<int, std::shared_ptr<PatchInfo<D>>>;
	/**
	 * @brief Use Zoltan to distribute patches
	 *
	 * The is only used on the coarsest domain
	 *
	 * @param level the level
	 */
	void balanceLevel(PInfoMap &level);
	/**
	 * @brief Balance a level with zoltan using a finer level
	 *
	 * @param level the new coarse level
	 * @param lower_level the finer level
	 */
	void balanceLevelWithLower(PInfoMap &level, PInfoMap &lower_level);

	public:
	/**
	 * @brief Construct a new DomGen object
	 *
	 * @param t the tree to use
	 * @param ns the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch.
	 * @param neumann whether to use neumann boundary conditions
	 */
	DomGen(Tree<D> t, std::array<int, D> ns, int num_ghost_cells, bool neumann = false);
	std::shared_ptr<Domain<D>> getFinestDomain();
	bool                       hasCoarserDomain();
	std::shared_ptr<Domain<D>> getCoarserDomain();
};
template <size_t D>
DomGen<D>::DomGen(Tree<D> t, std::array<int, D> ns, int num_ghost_cells, bool neumann)
{
	this->t               = t;
	this->ns              = ns;
	this->neumann         = neumann;
	this->num_ghost_cells = num_ghost_cells;
	num_levels            = t.num_levels;
	curr_level            = num_levels;

	// generate finest DC
	extractLevel();
}
template <size_t D> std::shared_ptr<Domain<D>> DomGen<D>::getFinestDomain()
{
	return domain_list.front();
}
template <size_t D> std::shared_ptr<Domain<D>> DomGen<D>::getCoarserDomain()
{
	if (curr_level > 0) {
		extractLevel();
		return domain_list.back();
	} else {
		return nullptr;
	}
}
template <size_t D> bool DomGen<D>::hasCoarserDomain()
{
	return curr_level > 0;
}
template <size_t D> inline void DomGen<D>::extractLevel()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	PInfoMap new_level;
	if (rank == 0) {
		Node<D>         child = *t.levels.at(curr_level);
		std::deque<int> q;
		std::set<int>   qed;
		q.push_back(child.id);
		qed.insert(child.id);

		while (!q.empty()) {
			std::shared_ptr<PatchInfo<D>> d_ptr(new PatchInfo<D>());
			PatchInfo<D> &                pinfo = *d_ptr;
			pinfo.num_ghost_cells               = num_ghost_cells;
			Node<D> n                           = t.nodes.at(q.front());
			q.pop_front();

			pinfo.ns = ns;
			pinfo.id = n.id;
			for (size_t i = 0; i < D; i++) {
				pinfo.spacings[i] = n.lengths[i] / ns[i];
			}
			pinfo.starts       = n.starts;
			pinfo.refine_level = n.level;
			if (n.level < curr_level) {
				// will still exist at next coarsening
				pinfo.parent_id = n.id;
				if (curr_level != num_levels) { pinfo.child_ids[0] = n.id; }
			} else {
				pinfo.parent_id = n.parent;
				if (pinfo.parent_id != -1) {
					char orth_on_parent = 0;
					while (t.nodes.at(n.parent).child_id[orth_on_parent] != n.id) {
						orth_on_parent++;
					}
					pinfo.orth_on_parent = orth_on_parent;
				}
				if (!n.hasChildren() && curr_level != num_levels) {
					pinfo.child_ids[0] = n.id;
				} else if (n.hasChildren()) {
					pinfo.child_ids = n.child_id;
				}
			}
			// set child ranks
			if (pinfo.child_ids[0] != -1) {
				for (size_t i = 0; i < pinfo.child_ids.size(); i++) {
					if (pinfo.child_ids[i] != -1) {
						pinfo.child_ranks[i] = id_rank_map.at(pinfo.child_ids[i]);
					}
				}
			}

			// set and enqueue nbrs
			for (Side<D> s : Side<D>::getValues()) {
				if (n.nbrId(s) == -1 && n.parent != -1 && t.nodes.at(n.parent).nbrId(s) != -1) {
					Node<D> parent = t.nodes.at(n.parent);
					Node<D> nbr    = t.nodes.at(parent.nbrId(s));
					auto    octs   = Orthant<D>::getValuesOnSide(s);
					int     quad   = 0;
					while (parent.childId(octs[quad]) != n.id) {
						quad++;
					}
					pinfo.nbr_info[s.toInt()].reset(new CoarseNbrInfo<D>(nbr.id, quad));
					if (!qed.count(nbr.id)) {
						q.push_back(nbr.id);
						qed.insert(nbr.id);
					}
				} else if (n.level < curr_level && n.nbrId(s) != -1
				           && t.nodes.at(n.nbrId(s)).hasChildren()) {
					Node<D> nbr  = t.nodes.at(n.nbrId(s));
					auto    octs = Orthant<D>::getValuesOnSide(s.opposite());
					std::array<int, Orthant<D>::num_orthants / 2> nbr_ids;
					for (size_t i = 0; i < Orthant<D>::num_orthants / 2; i++) {
						int id     = nbr.childId(octs[i]);
						nbr_ids[i] = id;
						if (!qed.count(id)) {
							q.push_back(id);
							qed.insert(id);
						}
					}
					pinfo.nbr_info[s.toInt()].reset(new FineNbrInfo<D>(nbr_ids));
				} else if (n.nbrId(s) != -1) {
					int id = n.nbrId(s);
					if (!qed.count(id)) {
						q.push_back(id);
						qed.insert(id);
					}
					pinfo.nbr_info[s.toInt()].reset(new NormalNbrInfo<D>(id));
				}
			}
			new_level[pinfo.id] = d_ptr;
		}
		for (auto &p : new_level) {
			p.second->setPtrs(new_level);
		}
	}
	// balance
	if (curr_level == num_levels) {
		balanceLevel(new_level);
	} else {
		auto &old_level = domain_list.back()->getPatchInfoMap();
		balanceLevelWithLower(new_level, old_level);
		// update parent ranks
		std::map<int, int> id_rank_map;
		// get outgoing information
		std::map<int, std::set<int>> out_info;
		for (auto pair : new_level) {
			id_rank_map[pair.first] = rank;
			auto &pinfo             = pair.second;
			for (int i = 0; i < Orthant<D>::num_orthants; i++) {
				if (pinfo->child_ranks[i] != -1 && pinfo->child_ranks[i] != rank) {
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
	domain_list.push_back(std::shared_ptr<Domain<D>>(new Domain<D>(new_level)));
	if (neumann) {
		IsNeumannFunc<D> inf = [](Side<D>, const std::array<double, D> &,
		                          const std::array<double, D> &) { return true; };
		domain_list.back()->setNeumann(inf);
	}
	curr_level--;
}
template <size_t D> inline void DomGen<D>::balanceLevel(PInfoMap &level)
{
	struct Zoltan_Struct *zz = Zoltan_Create(MPI_COMM_WORLD);

	// parameters
	Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");  /* Zoltan method: HYPERGRAPH */
	Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");     /* don't use local IDs */
	Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");      /* we omit object weights */
	Zoltan_Set_Param(zz, "AUTO_MIGRATE", "FALSE");    /* we omit object weights */
	Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");         /* we omit object weights */
	Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");    /* we omit object weights */

	// Query functions
	// Number of Vertices
	auto numObjFn = [](void *data, int *ierr) -> int {
		PInfoMap &map = *(PInfoMap *) data;
		*ierr         = ZOLTAN_OK;
		return map.size();
	};
	Zoltan_Set_Num_Obj_Fn(zz, numObjFn, &level);

	// List of vertices
	auto objListFn
	= [](void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids,
	     ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr) {
		  PInfoMap &map = *(PInfoMap *) data;
		  *ierr         = ZOLTAN_OK;
		  int pos       = 0;
		  for (auto p : map) {
			  global_ids[pos] = p.first;
			  pos++;
		  }
	  };
	Zoltan_Set_Obj_List_Fn(zz, objListFn, &level);

	// Construct hypergraph
	struct CompressedVertex {
		std::vector<int> vertices;
		std::vector<int> ptrs;
		std::vector<int> edges;
	};
	CompressedVertex graph;
	for (auto &p : level) {
		PatchInfo<D> &pinfo = *p.second;
		graph.vertices.push_back(pinfo.id);
		graph.ptrs.push_back(graph.edges.size());
		for (Side<D> s : Side<D>::getValues()) {
			int edge_id = -1;
			if (pinfo.hasNbr(s)) {
				switch (pinfo.getNbrType(s)) {
					case NbrType::Normal:
						if (s.isLowerOnAxis()) {
							edge_id = pinfo.getNormalNbrInfo(s).id;
						} else {
							edge_id = pinfo.id;
						}
						break;
					case NbrType::Fine:
						edge_id = pinfo.id;
						break;
					case NbrType::Coarse:
						edge_id = pinfo.getCoarseNbrInfo(s).id;
						break;
				}
			}
			graph.edges.push_back(edge_id * Side<D>::num_sides + s.toInt());
		}
	}

	// set graph functions
	Zoltan_Set_HG_Size_CS_Fn(zz,
	                         [](void *data, int *num_lists, int *num_pins, int *format, int *ierr) {
		                         CompressedVertex &graph = *(CompressedVertex *) data;
		                         *ierr                   = ZOLTAN_OK;
		                         *num_lists              = graph.vertices.size();
		                         *num_pins               = graph.edges.size();
		                         *format                 = ZOLTAN_COMPRESSED_VERTEX;
	                         },
	                         &graph);
	Zoltan_Set_HG_CS_Fn(zz,
	                    [](void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
	                       int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
	                       ZOLTAN_ID_PTR pin_GID, int *ierr) {
		                    CompressedVertex &graph = *(CompressedVertex *) data;
		                    *ierr                   = ZOLTAN_OK;
		                    copy(graph.vertices.begin(), graph.vertices.end(), vtxedge_GID);
		                    copy(graph.ptrs.begin(), graph.ptrs.end(), vtxedge_ptr);
		                    copy(graph.edges.begin(), graph.edges.end(), pin_GID);
	                    },
	                    &graph);

	Zoltan_Set_Obj_Size_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr) {
		                       PInfoMap &map = *(PInfoMap *) data;
		                       *ierr         = ZOLTAN_OK;
		                       return map[*global_id]->serialize(nullptr);
	                       },
	                       &level);
	Zoltan_Set_Pack_Obj_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size,
	                          char *buf, int *ierr) {
		                       PInfoMap &map = *(PInfoMap *) data;
		                       *ierr         = ZOLTAN_OK;
		                       map[*global_id]->serialize(buf);
		                       map.erase(*global_id);
	                       },
	                       &level);
	Zoltan_Set_Unpack_Obj_Fn(
	zz,
	[](void *data, int num_gid_entries, ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr) {
		PInfoMap &map = *(PInfoMap *) data;
		*ierr         = ZOLTAN_OK;
		map[*global_id].reset(new PatchInfo<D>());
		map[*global_id]->deserialize(buf);
	},
	&level);
	// Zoltan_Set_Obj_Size_Multi_Fn(zz, DomainZoltanHelper::object_sizes, this);
	// Zoltan_Set_Pack_Obj_Multi_Fn(zz, DomainZoltanHelper::pack_objects, this);
	// Zoltan_Set_Unpack_Obj_Multi_Fn(zz, DomainZoltanHelper::unpack_objects, this);

	////////////////////////////////////////////////////////////////
	// Zoltan can now partition the objects in this collection.
	// In this simple example, we assume the number of partitions is
	// equal to the number of processes.  Process rank 0 will own
	// partition 0, process rank 1 will own partition 1, and so on.
	////////////////////////////////////////////////////////////////

	int           changes;
	int           numGidEntries;
	int           numLidEntries;
	int           numImport;
	ZOLTAN_ID_PTR importGlobalIds;
	ZOLTAN_ID_PTR importLocalIds;
	int *         importProcs;
	int *         importToPart;
	int           numExport;
	ZOLTAN_ID_PTR exportGlobalIds;
	ZOLTAN_ID_PTR exportLocalIds;
	int *         exportProcs;
	int *         exportToPart;

	int rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries, &numImport,
	                             &importGlobalIds, &importLocalIds, &importProcs, &importToPart,
	                             &numExport, &exportGlobalIds, &exportLocalIds, &exportProcs,
	                             &exportToPart);

	// update ranks of neighbors before migrating
	id_rank_map.clear();
	for (int i = 0; i < numExport; i++) {
		int curr_id          = exportGlobalIds[i];
		int dest_rank        = exportProcs[i];
		id_rank_map[curr_id] = dest_rank;
		level[curr_id]->updateRank(dest_rank);
	}

	rc = Zoltan_Migrate(zz, numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
	                    numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		std::cerr << "zoltan error\n";
		Zoltan_Destroy(&zz);
		exit(0);
	}
	Zoltan_Destroy(&zz);
#if DD_DEBUG
	std::cout << "I have " << level.size() << " domains: ";

	int  prev  = -100;
	bool range = false;
	for (auto &p : level) {
		int curr = p.second->id;
		if (curr != prev + 1 && !range) {
			std::cout << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			std::cout << prev << " " << curr << "-";
		}
		prev = curr;
	}

	std::cout << prev << "\n";
	std::cout << std::endl;
#endif
	for (auto &p : level) {
		p.second->setPtrs(level);
	}
}
template <size_t D>
inline void DomGen<D>::balanceLevelWithLower(PInfoMap &level, PInfoMap &lower_level)
{
	struct Zoltan_Struct *zz = Zoltan_Create(MPI_COMM_WORLD);

	// parameters
	Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");  /* Zoltan method: HYPERGRAPH */
	Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");     /* don't use local IDs */
	Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");      /* we omit object weights */
	Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");     /* we omit object weights */
	Zoltan_Set_Param(zz, "AUTO_MIGRATE", "FALSE");    /* we omit object weights */
	Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");         /* we omit object weights */
	Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");    /* we omit object weights */

	// Query functions
	// Number of Vertices
	struct Levels {
		PInfoMap *upper;
		PInfoMap *lower;
	};
	Levels levels = {&level, &lower_level};
	Zoltan_Set_Num_Obj_Fn(zz,
	                      [](void *data, int *ierr) -> int {
		                      Levels *levels = (Levels *) data;
		                      *ierr          = ZOLTAN_OK;
		                      return levels->upper->size() + levels->lower->size();
	                      },
	                      &levels);

	// List of vertices
	Zoltan_Set_Obj_List_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim,
	                          float *obj_wgts, int *ierr) {
		                       Levels *levels = (Levels *) data;
		                       *ierr          = ZOLTAN_OK;
		                       int pos        = 0;
		                       for (auto p : *levels->upper) {
			                       global_ids[pos] = p.first;
			                       obj_wgts[pos]   = 1;
			                       pos++;
		                       }
		                       for (auto p : *levels->lower) {
			                       global_ids[pos] = p.first + levels->upper->size();
			                       obj_wgts[pos]   = 0;
			                       pos++;
		                       }
	                       },
	                       &levels);

	// Construct hypergraph
	struct CompressedVertex {
		std::vector<int> vertices;
		std::vector<int> ptrs;
		std::vector<int> edges;
	};
	CompressedVertex graph;
	// process coarse level
	for (auto &p : level) {
		PatchInfo<D> &pinfo = *p.second;
		graph.vertices.push_back(pinfo.id);
		graph.ptrs.push_back(graph.edges.size());
		// patch to patch communication
		for (Side<D> s : Side<D>::getValues()) {
			int edge_id = -1;
			if (pinfo.hasNbr(s)) {
				switch (pinfo.getNbrType(s)) {
					case NbrType::Normal:
						if (s.isLowerOnAxis()) {
							edge_id = pinfo.getNormalNbrInfo(s).id;
						} else {
							edge_id = pinfo.id;
						}
						break;
					case NbrType::Fine:
						edge_id = pinfo.id;
						break;
					case NbrType::Coarse:
						edge_id = pinfo.getCoarseNbrInfo(s).id;
						break;
				}
			}
			graph.edges.push_back(edge_id * Side<D>::num_sides + s.toInt());
		}
		// level to level communication
		graph.edges.push_back(-pinfo.id - 1);
	}
	// process fine level
	for (auto &p : lower_level) {
		PatchInfo<D> &pinfo = *p.second;
		graph.vertices.push_back(pinfo.id + level.size());
		graph.ptrs.push_back(graph.edges.size());
		graph.edges.push_back(-pinfo.parent_id - 1);
	}

	// set graph functions
	Zoltan_Set_HG_Size_CS_Fn(zz,
	                         [](void *data, int *num_lists, int *num_pins, int *format, int *ierr) {
		                         CompressedVertex &graph = *(CompressedVertex *) data;
		                         *ierr                   = ZOLTAN_OK;
		                         *num_lists              = graph.vertices.size();
		                         *num_pins               = graph.edges.size();
		                         *format                 = ZOLTAN_COMPRESSED_VERTEX;
	                         },
	                         &graph);
	Zoltan_Set_HG_CS_Fn(zz,
	                    [](void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
	                       int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
	                       ZOLTAN_ID_PTR pin_GID, int *ierr) {
		                    CompressedVertex &graph = *(CompressedVertex *) data;
		                    *ierr                   = ZOLTAN_OK;
		                    copy(graph.vertices.begin(), graph.vertices.end(), vtxedge_GID);
		                    copy(graph.ptrs.begin(), graph.ptrs.end(), vtxedge_ptr);
		                    copy(graph.edges.begin(), graph.edges.end(), pin_GID);
	                    },
	                    &graph);

	Zoltan_Set_Obj_Size_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr) {
		                       Levels *levels = (Levels *) data;
		                       *ierr          = ZOLTAN_OK;
		                       return levels->upper->at(*global_id)->serialize(nullptr);
	                       },
	                       &levels);

	// fixed objects
	Zoltan_Set_Num_Fixed_Obj_Fn(zz,
	                            [](void *data, int *ierr) -> int {
		                            Levels *levels = (Levels *) data;
		                            *ierr          = ZOLTAN_OK;
		                            return levels->lower->size();
	                            },
	                            &levels);

	Zoltan_Set_Fixed_Obj_List_Fn(zz,
	                             [](void *data, int num_fixed_obj, int num_gid_entries,
	                                ZOLTAN_ID_PTR fixed_gids, int *fixed_parts, int *ierr) {
		                             Levels *levels = (Levels *) data;
		                             *ierr          = ZOLTAN_OK;
		                             int rank;
		                             MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		                             int pos = 0;
		                             for (auto p : *levels->lower) {
			                             fixed_gids[pos]  = p.first + levels->upper->size();
			                             fixed_parts[pos] = rank;
			                             pos++;
		                             }
	                             },
	                             &levels);
	// pack and unpack
	// pack and unpack
	Zoltan_Set_Pack_Obj_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size,
	                          char *buf, int *ierr) {
		                       Levels *levels = (Levels *) data;
		                       *ierr          = ZOLTAN_OK;
		                       levels->upper->at(*global_id)->serialize(buf);
		                       levels->upper->erase(*global_id);
	                       },
	                       &levels);
	Zoltan_Set_Unpack_Obj_Fn(
	zz,
	[](void *data, int num_gid_entries, ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr) {
		Levels *levels = (Levels *) data;
		*ierr          = ZOLTAN_OK;
		(*levels->upper)[*global_id].reset(new PatchInfo<D>());
		(*levels->upper)[*global_id]->deserialize(buf);
	},
	&levels);

	////////////////////////////////////////////////////////////////
	// Zoltan can now partition the objects in this collection.
	// In this simple example, we assume the number of partitions is
	// equal to the number of processes.  Process rank 0 will own
	// partition 0, process rank 1 will own partition 1, and so on.
	////////////////////////////////////////////////////////////////

	int           changes;
	int           numGidEntries;
	int           numLidEntries;
	int           numImport;
	ZOLTAN_ID_PTR importGlobalIds;
	ZOLTAN_ID_PTR importLocalIds;
	int *         importProcs;
	int *         importToPart;
	int           numExport;
	ZOLTAN_ID_PTR exportGlobalIds;
	ZOLTAN_ID_PTR exportLocalIds;
	int *         exportProcs;
	int *         exportToPart;

	int rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries, &numImport,
	                             &importGlobalIds, &importLocalIds, &importProcs, &importToPart,
	                             &numExport, &exportGlobalIds, &exportLocalIds, &exportProcs,
	                             &exportToPart);

	id_rank_map.clear();
	// update ranks of neighbors before migrating
	for (int i = 0; i < numExport; i++) {
		int  curr_id   = exportGlobalIds[i];
		int  dest_rank = exportProcs[i];
		auto iter      = levels.upper->find(curr_id);
		if (iter != levels.upper->end()) {
			id_rank_map[curr_id] = dest_rank;
			iter->second->updateRank(dest_rank);
		}
	}

	rc = Zoltan_Migrate(zz, numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
	                    numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		std::cerr << "zoltan error\n";
		Zoltan_Destroy(&zz);
		exit(0);
	}
	Zoltan_Destroy(&zz);
#if DD_DEBUG
	std::cout << "I have " << levels.upper->size() << " domains: ";

	int  prev  = -100;
	bool range = false;
	for (auto &p : *levels.upper) {
		int curr = p.second->id;
		if (curr != prev + 1 && !range) {
			std::cout << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			std::cout << prev << " " << curr << "-";
		}
		prev = curr;
	}

	std::cout << prev << "\n";
	std::cout << std::endl;
#endif
	for (auto &p : level) {
		p.second->setPtrs(level);
	}
}
extern template class DomGen<2>;
extern template class DomGen<3>;
} // namespace Experimental
} // namespace Thunderegg
#endif
