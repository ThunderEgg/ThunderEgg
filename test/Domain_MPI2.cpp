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
#include <ThunderEgg/Domain.h>

#include "utils/DomainReader.h"

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Domain<3> numLocalPatches")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    CHECK_EQ(domain.getNumLocalPatches(), 14);
  } else {
    CHECK_EQ(domain.getNumLocalPatches(), 1);
  }
}
TEST_CASE("Domain<3> numGlobalPatches")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  CHECK_EQ(domain.getNumGlobalPatches(), 15);
}
TEST_CASE("Domain<3> numLocalCells")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    CHECK_EQ(domain.getNumLocalCells(), 14 * 10 * 10 * 10);
  } else {
    CHECK_EQ(domain.getNumLocalCells(), 1 * 10 * 10 * 10);
  }
}
TEST_CASE("Domain<3> numGlobalCells")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  CHECK_EQ(domain.getNumGlobalCells(), 15 * 10 * 10 * 10);
}
TEST_CASE("Domain<3> getNumGhostCells")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 2);
  auto domain = domain_reader.getFinerDomain();
  CHECK_EQ(domain.getNumGhostCells(), 2);
}
TEST_CASE("Domain<3> numLocalCellsWithGhost")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 2);
  auto domain = domain_reader.getFinerDomain();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    CHECK_EQ(domain.getNumLocalCellsWithGhost(), 14 * 14 * 14 * 14);
  } else {
    CHECK_EQ(domain.getNumLocalCellsWithGhost(), 1 * 14 * 14 * 14);
  }
}
TEST_CASE("Domain<3> getNs")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 11, 12 }, 2);
  auto domain = domain_reader.getFinerDomain();

  CHECK_EQ(domain.getNs()[0], 10);
  CHECK_EQ(domain.getNs()[1], 11);
  CHECK_EQ(domain.getNs()[2], 12);
}
TEST_CASE("Domain<3> getPatchInfoVector size")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    CHECK_EQ(domain.getPatchInfoVector().size(), 14);
  } else {
    CHECK_EQ(domain.getPatchInfoVector().size(), 1);
  }
}
TEST_CASE("Domain<3> local indexes match position in pinfo vector")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  auto pinfo_vector = domain.getPatchInfoVector();
  for (size_t i = 0; i < pinfo_vector.size(); i++) {
    CHECK_EQ(pinfo_vector[i].local_index, (int)i);
  }
}
TEST_CASE("Schur::InterfaceDomain<3> local indexes in neighbor info are consistent")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  map<int, int> id_to_local_index_map;

  for (auto pinfo : domain.getPatchInfoVector()) {
    id_to_local_index_map[pinfo.id] = pinfo.local_index;
  }

  auto checkIdAndLocalIndex = [&](int id, int local_index) {
    auto iter = id_to_local_index_map.find(id);
    if (iter != id_to_local_index_map.end()) {
      CHECK_EQ(local_index, iter->second);
    } else {
      CHECK_EQ(local_index, -1);
    }
  };
  for (auto pinfo : domain.getPatchInfoVector()) {
    for (Side<3> s : Side<3>::getValues()) {
      if (pinfo.hasNbr(s)) {
        NbrType type = pinfo.getNbrType(s);
        if (type == NbrType::Normal) {
          const NormalNbrInfo<2>& info = pinfo.getNormalNbrInfo(s);
          checkIdAndLocalIndex(info.id, info.local_index);
        } else if (type == NbrType::Fine) {
          const FineNbrInfo<2>& info = pinfo.getFineNbrInfo(s);
          for (size_t i = 0; i < info.ids.size(); i++) {
            checkIdAndLocalIndex(info.ids[i], info.local_indexes[i]);
          }
        } else if (type == NbrType::Coarse) {
          const CoarseNbrInfo<2>& info = pinfo.getCoarseNbrInfo(s);
          checkIdAndLocalIndex(info.id, info.local_index);
        }
      }
    }
    for (Edge e : Edge::getValues()) {
      if (pinfo.hasNbr(e)) {
        NbrType type = pinfo.getNbrType(e);
        if (type == NbrType::Normal) {
          const NormalNbrInfo<1>& info = pinfo.getNormalNbrInfo(e);
          checkIdAndLocalIndex(info.id, info.local_index);
        } else if (type == NbrType::Fine) {
          const FineNbrInfo<1>& info = pinfo.getFineNbrInfo(e);
          for (size_t i = 0; i < info.ids.size(); i++) {
            checkIdAndLocalIndex(info.ids[i], info.local_indexes[i]);
          }
        } else if (type == NbrType::Coarse) {
          const CoarseNbrInfo<1>& info = pinfo.getCoarseNbrInfo(e);
          checkIdAndLocalIndex(info.id, info.local_index);
        }
      }
    }
    for (Corner<3> c : Corner<3>::getValues()) {
      if (pinfo.hasNbr(c)) {
        NbrType type = pinfo.getNbrType(c);
        if (type == NbrType::Normal) {
          const NormalNbrInfo<0>& info = pinfo.getNormalNbrInfo(c);
          checkIdAndLocalIndex(info.id, info.local_index);
        } else if (type == NbrType::Fine) {
          const FineNbrInfo<0>& info = pinfo.getFineNbrInfo(c);
          for (size_t i = 0; i < info.ids.size(); i++) {
            checkIdAndLocalIndex(info.ids[i], info.local_indexes[i]);
          }
        } else if (type == NbrType::Coarse) {
          const CoarseNbrInfo<0>& info = pinfo.getCoarseNbrInfo(c);
          checkIdAndLocalIndex(info.id, info.local_index);
        }
      }
    }
  }
}
TEST_CASE("Domain<3> global indexes match position in pinfo vector")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int start_i;
  if (rank == 0) {
    start_i = 0;
  } else {
    start_i = 14;
  }
  auto pinfo_vector = domain.getPatchInfoVector();
  for (size_t i = 0; i < pinfo_vector.size(); i++) {
    CHECK_EQ(pinfo_vector[i].global_index, start_i + (int)i);
  }
}
TEST_CASE("Domain<3> global indexes in neighbor info are consistent")
{
  DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json", { 10, 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();

  std::vector<int> global_index_to_id_in(domain.getNumGlobalPatches());
  std::vector<int> global_index_to_id_out(domain.getNumGlobalPatches());

  for (auto pinfo : domain.getPatchInfoVector()) {
    global_index_to_id_in[pinfo.global_index] = pinfo.id;
  }

  MPI_Allreduce(global_index_to_id_in.data(), global_index_to_id_out.data(), global_index_to_id_in.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  map<int, int> id_to_global_index_map;
  for (size_t i = 0; i < global_index_to_id_out.size(); i++) {
    id_to_global_index_map[global_index_to_id_out[i]] = i;
  }

  auto checkIdAndLocalIndex = [&](int id, int global_index) { CHECK_EQ(global_index, id_to_global_index_map.at(id)); };
  for (auto pinfo : domain.getPatchInfoVector()) {
    for (Side<3> s : Side<3>::getValues()) {
      if (pinfo.hasNbr(s)) {
        NbrType type = pinfo.getNbrType(s);
        if (type == NbrType::Normal) {
          const NormalNbrInfo<2>& info = pinfo.getNormalNbrInfo(s);
          checkIdAndLocalIndex(info.id, info.global_index);
        } else if (type == NbrType::Fine) {
          const FineNbrInfo<2>& info = pinfo.getFineNbrInfo(s);
          for (size_t i = 0; i < info.ids.size(); i++) {
            checkIdAndLocalIndex(info.ids[i], info.global_indexes[i]);
          }
        } else if (type == NbrType::Coarse) {
          const CoarseNbrInfo<2>& info = pinfo.getCoarseNbrInfo(s);
          checkIdAndLocalIndex(info.id, info.global_index);
        }
      }
    }
    for (Edge e : Edge::getValues()) {
      if (pinfo.hasNbr(e)) {
        NbrType type = pinfo.getNbrType(e);
        if (type == NbrType::Normal) {
          const NormalNbrInfo<1>& info = pinfo.getNormalNbrInfo(e);
          checkIdAndLocalIndex(info.id, info.global_index);
        } else if (type == NbrType::Fine) {
          const FineNbrInfo<1>& info = pinfo.getFineNbrInfo(e);
          for (size_t i = 0; i < info.ids.size(); i++) {
            checkIdAndLocalIndex(info.ids[i], info.global_indexes[i]);
          }
        } else if (type == NbrType::Coarse) {
          const CoarseNbrInfo<1>& info = pinfo.getCoarseNbrInfo(e);
          checkIdAndLocalIndex(info.id, info.global_index);
        }
      }
    }
    for (Corner<3> c : Corner<3>::getValues()) {
      if (pinfo.hasNbr(c)) {
        NbrType type = pinfo.getNbrType(c);
        if (type == NbrType::Normal) {
          const NormalNbrInfo<0>& info = pinfo.getNormalNbrInfo(c);
          checkIdAndLocalIndex(info.id, info.global_index);
        } else if (type == NbrType::Fine) {
          const FineNbrInfo<0>& info = pinfo.getFineNbrInfo(c);
          for (size_t i = 0; i < info.ids.size(); i++) {
            checkIdAndLocalIndex(info.ids[i], info.global_indexes[i]);
          }
        } else if (type == NbrType::Coarse) {
          const CoarseNbrInfo<0>& info = pinfo.getCoarseNbrInfo(c);
          checkIdAndLocalIndex(info.id, info.global_index);
        }
      }
    }
  }
}
