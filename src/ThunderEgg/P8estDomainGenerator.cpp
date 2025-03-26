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

#include "P8estDomainGenerator.h"

#include <ThunderEgg/CoarseNbrInfo.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/Face.h>
#include <ThunderEgg/FineNbrInfo.h>
#include <ThunderEgg/NormalNbrInfo.h>
#include <ThunderEgg/Orthant.h>
#include <ThunderEgg/PatchInfo.h>
#include <ThunderEgg/RuntimeError.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <mpi.h>
#include <p4est_base.h>
#include <p8est.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_iterate.h>
#include <p8est_mesh.h>
#include <sc.h>
#include <sc_containers.h>
#include <set>
#include <vector>

using namespace std;
using namespace ThunderEgg;

namespace {
struct IterateFunctions
{
  std::function<void(p8est_iter_volume_info_t*)> iter_volume;
  std::function<void(p8est_iter_face_info_t*)> iter_face;
  std::function<void(p8est_iter_edge_info_t*)> iter_edge;
  std::function<void(p8est_iter_corner_info_t*)> iter_corner;
};

void
IterVolumeWrap(p8est_iter_volume_info_t* info, void* user_data)
{
  const IterateFunctions* functions = (IterateFunctions*)user_data;
  functions->iter_volume(info);
}

void
IterFaceWrap(p8est_iter_face_info_t* info, void* user_data)
{
  const IterateFunctions* functions = (IterateFunctions*)user_data;
  functions->iter_face(info);
}

void
IterEdgeWrap(p8est_iter_edge_info_t* info, void* user_data)
{
  const IterateFunctions* functions = (IterateFunctions*)user_data;
  functions->iter_edge(info);
}

void
IterCornerWrap(p8est_iter_corner_info_t* info, void* user_data)
{
  const IterateFunctions* functions = (IterateFunctions*)user_data;
  functions->iter_corner(info);
}

void
p8est_iterate_ext(p8est_t* p8est, p8est_ghost_t* ghost_layer, const std::function<void(const p8est_iter_volume_info_t*)>& iter_volume, const std::function<void(const p8est_iter_face_info_t*)>& iter_face, const std::function<void(const p8est_iter_edge_info_t*)>& iter_edge, const std::function<void(const p8est_iter_corner_info_t*)>& iter_corner, int remote)
{
  IterateFunctions functions = { iter_volume, iter_face, iter_edge, iter_corner };
  p8est_iterate_ext(p8est, ghost_layer, &functions, iter_volume ? &IterVolumeWrap : nullptr, iter_face ? &IterFaceWrap : nullptr, iter_edge ? &IterEdgeWrap : nullptr, iter_corner ? &IterCornerWrap : nullptr, remote);
}

struct CoarsenFunctions
{
  std::function<int(p4est_topidx_t, p8est_quadrant_t*[])> coarsen;
  std::function<void(p4est_topidx_t, int, p8est_quadrant_t*[], int, p8est_quadrant_t*[])> replace;
};

int
CoarsenWrap(p8est_t* p8est, p4est_topidx_t which_tree, p8est_quadrant_t* quadrants[])
{
  const CoarsenFunctions* functions = (CoarsenFunctions*)p8est->user_pointer;
  return functions->coarsen(which_tree, quadrants);
}

void
ReplaceWrap(p8est_t* p8est, p4est_topidx_t which_tree, int num_outgoing, p8est_quadrant_t* outgoing[], int num_incoming, p8est_quadrant_t* incoming[])
{
  const CoarsenFunctions* functions = (CoarsenFunctions*)p8est->user_pointer;
  functions->replace(which_tree, num_outgoing, outgoing, num_incoming, incoming);
}

void
p8est_coarsen_ext(p8est_t* p8est, int coarsen_recursive, int callback_orphans, const std::function<int(p4est_topidx_t, p8est_quadrant_t*[])>& coarsen, p8est_init_t init_fn, const std::function<void(p4est_topidx_t, int, p8est_quadrant_t*[], int, p8est_quadrant_t*[])>& replace)
{
  CoarsenFunctions functions = { coarsen, replace };
  p8est->user_pointer = &functions;
  p8est_coarsen_ext(p8est, coarsen_recursive, callback_orphans, coarsen ? &CoarsenWrap : nullptr, init_fn, replace ? &ReplaceWrap : nullptr);
  p8est->user_pointer = nullptr;
}

int
GetMaxLevel(p8est_t* p8est)
{
  int max_level = 0;

  p8est_iterate_ext(p8est, nullptr, [&](const p8est_iter_volume_info_t* info) { max_level = max(max_level, (int)info->quad->level); }, nullptr, nullptr, nullptr, false);

  int global_max_level;

  MPI_Allreduce(&max_level, &global_max_level, 1, MPI_INT, MPI_MAX, p8est->mpicomm);
  return global_max_level;
}

struct Data
{
  int id;
  int rank;
  int level;
  std::array<int, 8> child_ids;
  std::array<int, 8> child_ranks;
  PatchInfo<3>* pinfo;
};

void
InitData(p8est_t*, p4est_topidx_t, p8est_quadrant_t* quadrant)
{
  Data& data = *(Data*)quadrant->p.user_data;
  data.id = -1;
  data.rank = -1;
  data.level = -1;
  data.child_ids.fill(-1);
  data.child_ranks.fill(-1);
  data.pinfo = nullptr;
}

void
SetRanks(p8est_t* p8est)
{
  int rank;
  MPI_Comm_rank(p8est->mpicomm, &rank);

  p8est_iterate_ext(
    p8est,
    nullptr,
    [=](const p8est_iter_volume_info_t* info) {
      Data* data = (Data*)info->quad->p.user_data;
      data->rank = rank;
    },
    nullptr,
    nullptr,
    nullptr,
    false);
}

void
SetIds(p8est_t* p8est)
{
  int curr_global_id;
  MPI_Scan(&p8est->local_num_quadrants, &curr_global_id, 1, MPI_INT, MPI_SUM, p8est->mpicomm);
  curr_global_id -= p8est->local_num_quadrants;

  p8est_iterate_ext(
    p8est,
    nullptr,
    [&](const p8est_iter_volume_info_t* info) {
      Data* data = (Data*)info->quad->p.user_data;
      data->id = curr_global_id;
      curr_global_id++;
    },
    nullptr,
    nullptr,
    nullptr,
    false);
}

void
SetLevels(p8est_t* p8est)
{
  p8est_iterate_ext(
    p8est,
    nullptr,
    [](const p8est_iter_volume_info_t* info) {
      Data* data = (Data*)info->quad->p.user_data;
      data->level = info->quad->level;
    },
    nullptr,
    nullptr,
    nullptr,
    false);
}

} // namespace

/*
 * constructor
 */
P8estDomainGenerator::P8estDomainGenerator(p8est_t* p8est, const std::array<int, 3>& ns, int num_ghost_cells, const BlockMapFunc& bmf)
  : ns(ns)
  , num_ghost_cells(num_ghost_cells)
  , bmf(bmf)
  , comm(p8est->mpicomm)
{
  my_p8est = p8est_copy(p8est, false);

  p8est_reset_data(my_p8est, sizeof(Data), InitData, nullptr);

  curr_level = GetMaxLevel(my_p8est);

  SetIds(my_p8est);

  extractLevel();
}

P8estDomainGenerator::~P8estDomainGenerator()
{
  p8est_destroy(my_p8est);
}

P8estDomainGenerator::P8estDomainGenerator(P8estDomainGenerator& other)
  : domain_patches(other.domain_patches)
  , ns(other.ns)
  , num_ghost_cells(other.num_ghost_cells)
  , curr_level(other.curr_level)
  , bmf(other.bmf)
{
  my_p8est = p8est_copy(other.my_p8est, true);
}

P8estDomainGenerator&
P8estDomainGenerator::operator=(const P8estDomainGenerator& other)
{
  domain_patches = other.domain_patches;
  ns = other.ns;
  num_ghost_cells = other.num_ghost_cells;
  curr_level = other.curr_level;
  bmf = other.bmf;
  p8est_destroy(my_p8est);
  my_p8est = p8est_copy(other.my_p8est, true);
  return *this;
}
void
P8estDomainGenerator::extractLevel()
{
  if (domain_patches.size() > 0) {
    coarsenTree();
  }

  SetRanks(my_p8est);
  SetLevels(my_p8est);

  createPatchInfos();

  linkNeighbors();

  if (domain_patches.size() == 2) {
    updateParentRanksOfPreviousDomain();
  }

  curr_level--;
}
void
P8estDomainGenerator::coarsenTree()
{
  p8est_coarsen_ext(
    my_p8est,
    false,
    true,
    [&](p4est_topidx_t, p8est_quadrant_t* quadrant[]) {
      if (quadrant[1] == nullptr) {
        // update child data to point to self
        Data* data = (Data*)quadrant[0]->p.user_data;
        data->child_ids[0] = data->id;
        data->child_ranks[0] = data->rank;
        data->pinfo->parent_id = data->id;
        data->pinfo->orth_on_parent = Orthant<3>::null();
        return false;
      } else if (quadrant[0]->level > curr_level) {
        // update child datas to point to parent
        const Data* bsw_data = (Data*)quadrant[0]->p.user_data;
        for (unsigned char i = 0; i < 8; i++) {
          Data* data = (Data*)quadrant[i]->p.user_data;
          data->child_ids[0] = data->id;
          data->child_ranks[0] = data->rank;
          data->pinfo->parent_id = bsw_data->id;
          data->pinfo->orth_on_parent = Orthant<3>(i);
        }
        return true;
      } else {
        // update child datas to point to self
        for (unsigned char i = 0; i < 8; i++) {
          Data* data = (Data*)quadrant[i]->p.user_data;
          data->child_ids[0] = data->id;
          data->child_ranks[0] = data->rank;
          data->pinfo->parent_id = data->id;
          data->pinfo->orth_on_parent = Orthant<3>::null();
        }
        return false;
      }
    },
    InitData,
    [](p4est_topidx_t, int, p8est_quadrant_t* outgoing[], int, p8est_quadrant_t* incoming[]) {
      Data* data = (Data*)incoming[0]->p.user_data;
      const Data* fine_sw_data = (Data*)outgoing[0]->p.user_data;
      data->id = fine_sw_data->id;
      for (int i = 0; i < 8; i++) {
        const Data* fine_data = (Data*)outgoing[i]->p.user_data;
        data->child_ids[i] = fine_data->id;
        data->child_ranks[i] = fine_data->rank;
      }
    });
  p8est_partition_ext(my_p8est, true, nullptr);
}

void
P8estDomainGenerator::createPatchInfos()
{
  domain_patches.emplace_front(my_p8est->local_num_quadrants);
  p8est_iterate_ext(
    my_p8est,
    nullptr,
    [&](const p8est_iter_volume_info_t* info) {
      Data* data = (Data*)info->quad->p.user_data;
      data->pinfo = &domain_patches.front()[info->quadid + p8est_tree_array_index(info->p4est->trees, info->treeid)->quadrants_offset];

      data->pinfo->rank = info->p4est->mpirank;
      data->pinfo->id = data->id;
      data->pinfo->ns = ns;
      data->pinfo->num_ghost_cells = num_ghost_cells;
      data->pinfo->refine_level = info->quad->level;

      data->pinfo->child_ids = data->child_ids;
      data->pinfo->child_ranks = data->child_ranks;

      double unit_x = (double)info->quad->x / P8EST_ROOT_LEN;
      double unit_y = (double)info->quad->y / P8EST_ROOT_LEN;
      double unit_z = (double)info->quad->z / P8EST_ROOT_LEN;

      bmf(info->treeid, unit_x, unit_y, unit_z, data->pinfo->starts[0], data->pinfo->starts[1], data->pinfo->starts[2]);

      double upper_unit_x = (double)(info->quad->x + P8EST_QUADRANT_LEN(info->quad->level)) / P8EST_ROOT_LEN;
      double upper_unit_y = (double)(info->quad->y + P8EST_QUADRANT_LEN(info->quad->level)) / P8EST_ROOT_LEN;
      double upper_unit_z = (double)(info->quad->z + P8EST_QUADRANT_LEN(info->quad->level)) / P8EST_ROOT_LEN;

      bmf(info->treeid, upper_unit_x, upper_unit_y, upper_unit_z, data->pinfo->spacings[0], data->pinfo->spacings[1], data->pinfo->spacings[2]);

      for (int i = 0; i < 3; i++) {
        data->pinfo->spacings[i] -= data->pinfo->starts[i];
        data->pinfo->spacings[i] /= ns[i];
      }
    },
    nullptr,
    nullptr,
    nullptr,
    false);
}

/*
 * Functions for linkNeighbors
 */
namespace {
Edge
getEdge(int p8est_edge)
{
  Edge edge;
  switch (p8est_edge) {
    case 0:
      edge = Edge::bs();
      break;
    case 1:
      edge = Edge::bn();
      break;
    case 2:
      edge = Edge::ts();
      break;
    case 3:
      edge = Edge::tn();
      break;
    case 4:
      edge = Edge::bw();
      break;
    case 5:
      edge = Edge::be();
      break;
    case 6:
      edge = Edge::tw();
      break;
    case 7:
      edge = Edge::te();
      break;
    case 8:
      edge = Edge::sw();
      break;
    case 9:
      edge = Edge::se();
      break;
    case 10:
      edge = Edge::nw();
      break;
    case 11:
      edge = Edge::ne();
      break;
    default:
      edge = Edge::null();
  }
  return edge;
}

void
SetFaceNbrInfo(p8est_mesh_t* mesh, vector<Data*>& data_ptrs, p4est_locidx_t quadid, int s)
{
  Data* data = data_ptrs[quadid];
  Side<3> side(s);
  int index = quadid * 6 + s;
  p4est_locidx_t qtq = mesh->quad_to_quad[index];
  int8_t qtf = mesh->quad_to_face[index];

  if (qtq == quadid && qtf == s) {
    return;
  }

  if (qtf >= 0 && qtf <= 23) {

    NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(data_ptrs[qtq]->id);
    nbr_info->rank = data_ptrs[qtq]->rank;

    data->pinfo->setNbrInfo(side, nbr_info);

  } else if (qtf >= 24 && qtf <= 119) {

    CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>();
    nbr_info->id = data_ptrs[qtq]->id;
    nbr_info->rank = data_ptrs[qtq]->rank;

    int8_t nbr_subface = (qtf - 24) / 24;
    nbr_info->orth_on_coarse = Orthant<2>(nbr_subface);

    data->pinfo->setNbrInfo(side, nbr_info);

  } else if (qtf >= -24 && qtf <= -1) {

    FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>();
    const p4est_locidx_t* qth = (p4est_locidx_t*)sc_array_index(mesh->quad_to_half, qtq);

    for (int i = 0; i < 4; i++) {
      nbr_info->ids[i] = data_ptrs[qth[i]]->id;
      nbr_info->ranks[i] = data_ptrs[qth[i]]->rank;
    }

    data->pinfo->setNbrInfo(side, nbr_info);

  } else {
    throw RuntimeError("Invalid quad_to_face value");
  }
}

void
SetEdgeNbrInfo(p8est_mesh_t* mesh, vector<Data*>& data_ptrs, p4est_locidx_t quadid, int i)
{
  Data* data = data_ptrs[quadid];
  Edge edge = getEdge(i);
  int index = quadid * 12 + i;
  p4est_locidx_t qte = mesh->quad_to_edge[index];

  if (qte < 0) {
    return;
  }

  if (qte < mesh->local_num_quadrants + mesh->ghost_num_quadrants) {
    // inter-tree normal nbr

    NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(data_ptrs[qte]->id);
    nbr_info->rank = data_ptrs[qte]->rank;

    data->pinfo->setNbrInfo(edge, nbr_info);

  } else {
    p4est_locidx_t offset_index = qte - (mesh->local_num_quadrants + mesh->ghost_num_quadrants);
    p4est_locidx_t start_index = *(p4est_locidx_t*)sc_array_index(mesh->edge_offset, offset_index);
    p4est_locidx_t end_index = *(p4est_locidx_t*)sc_array_index(mesh->edge_offset, offset_index + 1);

    int8_t ee = *(int8_t*)sc_array_index(mesh->edge_edge, start_index);

    if ((ee >= 0 && end_index - start_index != 1) || (ee < 0 && end_index - start_index != 2)) {
      throw RuntimeError("Unsupported number of edge neighbors");
    }

    if (ee >= 0 && ee <= 23) {
      // intra-tree normal nbr

      p4est_locidx_t eq = *(p4est_locidx_t*)sc_array_index(mesh->edge_quad, start_index);

      NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(data_ptrs[eq]->id);
      nbr_info->rank = data_ptrs[eq]->rank;

      data->pinfo->setNbrInfo(edge, nbr_info);

    } else if (ee >= 24 && ee <= 71) {

      p4est_locidx_t eq = *(p4est_locidx_t*)sc_array_index(mesh->edge_quad, start_index);

      CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>();
      nbr_info->id = data_ptrs[eq]->id;
      nbr_info->rank = data_ptrs[eq]->rank;

      uint8_t nbr_subface = (ee - 24) / 24;
      nbr_info->orth_on_coarse = Orthant<1>(nbr_subface);

      data->pinfo->setNbrInfo(edge, nbr_info);

    } else if (ee >= -24 && ee <= -1) {

      FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>();

      for (int i = 0; i < 2; i++) {
        p4est_locidx_t eq = *(p4est_locidx_t*)sc_array_index(mesh->edge_quad, start_index + i);
        nbr_info->ids[i] = data_ptrs[eq]->id;
        nbr_info->ranks[i] = data_ptrs[eq]->rank;
      }

      data->pinfo->setNbrInfo(edge, nbr_info);

    } else {
      throw RuntimeError("Invalid quad_to_edge_face value");
    }
  }
}

void
SetCornerNbrInfo(p8est_mesh_t* mesh, vector<Data*>& data_ptrs, p4est_locidx_t quadid, int i)
{
  Data* data = data_ptrs[quadid];
  Corner<3> corner(i);
  int index = quadid * 8 + i;
  p4est_locidx_t qtc = mesh->quad_to_corner[index];

  if (qtc < 0) {
    return;
  }

  int nbr_id;
  int nbr_rank;
  int level_diff;
  if (qtc < mesh->local_num_quadrants + mesh->ghost_num_quadrants) {
    // inter-tree normal nbr

    nbr_id = data_ptrs[qtc]->id;
    nbr_rank = data_ptrs[qtc]->rank;
    level_diff = data->level - data_ptrs[qtc]->level;
  } else {
    p4est_locidx_t offset_index = qtc - (mesh->local_num_quadrants + mesh->ghost_num_quadrants);
    p4est_locidx_t start_index = *(p4est_locidx_t*)sc_array_index(mesh->corner_offset, offset_index);
    p4est_locidx_t end_index = *(p4est_locidx_t*)sc_array_index(mesh->corner_offset, offset_index + 1);

    if (end_index - start_index != 1) {
      throw RuntimeError("Unsupported number of corner neighbors");
    }

    p4est_locidx_t cq = *(p4est_locidx_t*)sc_array_index(mesh->corner_quad, start_index);
    nbr_id = data_ptrs[cq]->id;
    nbr_rank = data_ptrs[cq]->rank;
    level_diff = data->level - data_ptrs[cq]->level;
  }
  if (level_diff < -1 || level_diff > 1) {
    throw RuntimeError("Invalid level difference between corner and neighbor");
  }
  if (level_diff == 0) {
    // normal nbr
    NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(nbr_id);
    nbr_info->rank = nbr_rank;
    data->pinfo->setNbrInfo(corner, nbr_info);
  } else if (level_diff == 1) {
    // coarse nbr
    CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>();
    nbr_info->orth_on_coarse = Orthant<0>(0);
    nbr_info->id = nbr_id;
    nbr_info->rank = nbr_rank;
    data->pinfo->setNbrInfo(corner, nbr_info);
  } else if (level_diff == -1) {
    // fine nbr
    FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>();
    nbr_info->ids[0] = nbr_id;
    nbr_info->ranks[0] = nbr_rank;
    data->pinfo->setNbrInfo(corner, nbr_info);
  }
}

} // namespace

void
P8estDomainGenerator::linkNeighbors()
{
  p8est_ghost_t* ghost = p8est_ghost_new(my_p8est, P8EST_CONNECT_CORNER);

  vector<Data> ghost_data(ghost->ghosts.elem_count);
  p8est_ghost_exchange_data(my_p8est, ghost, ghost_data.data());

  p8est_mesh_params_t mesh_params;
  p8est_mesh_params_init(&mesh_params);
  mesh_params.edgehanging_corners = 1;
  mesh_params.btype = P8EST_CONNECT_CORNER;

  p8est_mesh_t* mesh = p8est_mesh_new_params(my_p8est, ghost, &mesh_params);

  vector<Data*> data_ptrs(my_p8est->local_num_quadrants + ghost->ghosts.elem_count);
  p4est_locidx_t curr_id = 0;
  p8est_iterate_ext(
    my_p8est,
    ghost,
    [&](const p8est_iter_volume_info_t* info) {
      data_ptrs[curr_id] = (Data*)info->quad->p.user_data;
      curr_id++;
    },
    nullptr,
    nullptr,
    nullptr,
    false);
  for (int i = 0; i < ghost->ghosts.elem_count; i++) {
    data_ptrs[my_p8est->local_num_quadrants + i] = &ghost_data[i];
  }

  for (p4est_locidx_t quadid = 0; quadid < mesh->local_num_quadrants; quadid++) {
    for (int i = 0; i < 6; i++) {
      SetFaceNbrInfo(mesh, data_ptrs, quadid, i);
    }

    for (int i = 0; i < 12; i++) {
      SetEdgeNbrInfo(mesh, data_ptrs, quadid, i);
    }

    for (int i = 0; i < 8; i++) {
      SetCornerNbrInfo(mesh, data_ptrs, quadid, i);
    }
  }

  p8est_ghost_destroy(ghost);
  p8est_mesh_destroy(mesh);
}

void
P8estDomainGenerator::updateParentRanksOfPreviousDomain()
{
  std::vector<PatchInfo<3>>& old_level = domain_patches.back();
  const std::vector<PatchInfo<3>>& new_level = domain_patches.front();

  std::set<int> local_new_level_ids;
  for (const auto& pinfo : new_level) {
    local_new_level_ids.insert(pinfo.id);
  }

  // update parent ranks
  std::map<int, int> id_rank_map;
  // get outgoing information
  std::map<int, std::set<int>> out_info;
  for (const PatchInfo<3>& pinfo : new_level) {
    id_rank_map[pinfo.id] = my_p8est->mpirank;
    for (int i = 0; i < 8; i++) {
      if (pinfo.child_ranks[i] != -1 && pinfo.child_ranks[i] != my_p8est->mpirank) {
        out_info[pinfo.child_ranks[i]].insert(pinfo.id);
      }
    }
  }
  // get incoming information
  std::set<int> incoming_ids;
  for (const PatchInfo<3>& pinfo : old_level) {
    if (!local_new_level_ids.count(pinfo.parent_id)) {
      incoming_ids.insert(pinfo.parent_id);
    }
  }
  // send info
  std::deque<std::vector<int>> buffers;
  std::vector<MPI_Request> send_requests;
  for (const auto& pair : out_info) {
    int dest = pair.first;
    std::vector<int> buffer(pair.second.begin(), pair.second.end());
    buffers.push_back(buffer);
    MPI_Request request;
    MPI_Isend(buffer.data(), (int)buffer.size(), MPI_INT, dest, 0, MPI_COMM_WORLD, &request);
    send_requests.push_back(request);
  }
  // recv info
  while (incoming_ids.size()) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    int size;
    MPI_Get_count(&status, MPI_INT, &size);
    int* buffer = new int[size];

    MPI_Recv(buffer, size, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);

    for (int i = 0; i < size; i++) {
      id_rank_map[buffer[i]] = status.MPI_SOURCE;
      incoming_ids.erase(buffer[i]);
    }

    delete[] buffer;
  }
  // wait for all
  MPI_Waitall((int)send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
  // update rank info
  for (PatchInfo<3>& pinfo : old_level) {
    pinfo.parent_rank = id_rank_map.at(pinfo.parent_id);
  }
}

Domain<3>
P8estDomainGenerator::getCoarserDomain()
{
  if (curr_level >= 0) {
    extractLevel();
  }
  Domain<3> domain(comm, id, ns, num_ghost_cells, domain_patches.back().begin(), domain_patches.back().end());
  domain_patches.pop_back();
  id++;
  return domain;
}

Domain<3>
P8estDomainGenerator::getFinestDomain()
{
  if (curr_level >= 0) {
    extractLevel();
  }
  Domain<3> domain(comm, id, ns, num_ghost_cells, domain_patches.back().begin(), domain_patches.back().end());
  domain_patches.pop_back();
  id++;
  return domain;
}

bool
P8estDomainGenerator::hasCoarserDomain()
{
  return !domain_patches.empty();
}
