/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#include "../utils/DomainReader.h"
#include <ThunderEgg/Schur/PatchIfaceScatter.h>

#include <limits>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

#define MESHES "mesh_inputs/2d_refined_complicated_mpi3.json"

/*****************
 *
 *
 *  Exception cases
 *
 *
 *****************/
TEST_CASE("Schur::PatchIfaceScatter<2> throws exception for non-square patches")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 5, 7 }) {
      for (auto ny : { 6, 8 }) {
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, 0);
        Domain<2> domain = domain_reader.getFinerDomain();
        Schur::InterfaceDomain<2> iface_domain(domain);

        auto construct = [&]() { Schur::PatchIfaceScatter<2> test(iface_domain); };
        CHECK_THROWS_AS(construct(), RuntimeError);
      }
    }
  }
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatterFinish throws exception when called with different "
          "global vectors")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      Vector<1> global_vector(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      Vector<1> global_vector_2(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      auto local_vector = scatter.getNewLocalPatchIfaceVector();

      auto state = scatter.scatterStart(global_vector, *local_vector);
      CHECK_THROWS_AS(scatter.scatterFinish(state, global_vector_2, *local_vector), RuntimeError);
    }
  }
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatterFinish throws exception when called with different "
          "local vectors")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      Vector<1> global_vector(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      auto local_vector = scatter.getNewLocalPatchIfaceVector();
      auto local_vector_2 = scatter.getNewLocalPatchIfaceVector();

      auto state = scatter.scatterStart(global_vector, *local_vector);
      CHECK_THROWS_AS(scatter.scatterFinish(state, global_vector, *local_vector_2), RuntimeError);
    }
  }
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatterFinish throws exception when called with different "
          "local and global vectors")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      Vector<1> global_vector(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      Vector<1> global_vector_2(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      auto local_vector = scatter.getNewLocalPatchIfaceVector();
      auto local_vector_2 = scatter.getNewLocalPatchIfaceVector();

      auto state = scatter.scatterStart(global_vector, *local_vector);
      CHECK_THROWS_AS(scatter.scatterFinish(state, global_vector_2, *local_vector_2), RuntimeError);
    }
  }
}
/******
 *
 *  getNewLocalPatchIfaceVector
 *
 *
 ********/
TEST_CASE("Schur::PatchIfaceScatter<2> getNewLocalPatchIfaceVector returns vector of expected length")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      auto local_vector = scatter.getNewLocalPatchIfaceVector();
      set<int> patch_iface_interfaces;

      for (auto piinfo : iface_domain.getPatchIfaceInfos()) {
        for (Side<2> s : Side<2>::getValues()) {
          if (piinfo->pinfo.hasNbr(s)) {
            auto iface_info = piinfo->getIfaceInfo(s);
            patch_iface_interfaces.insert(iface_info->id);
          }
        }
      }

      CHECK_EQ(local_vector->getNumLocalPatches(), (int)patch_iface_interfaces.size());
      CHECK_EQ(local_vector->getNumLocalCells(), n * (int)patch_iface_interfaces.size());
    }
  }
}
TEST_CASE("Schur::PatchIfaceScatter<2> getNewLocalPatchIfaceVector returns vector with local MPI_Comm")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      auto local_vector = scatter.getNewLocalPatchIfaceVector();

      int result;
      int err = MPI_Comm_compare(local_vector->getCommunicator().getMPIComm(), MPI_COMM_SELF, &result);
      REQUIRE_EQ(err, MPI_SUCCESS);
      CHECK_EQ(result, MPI_CONGRUENT);
    }
  }
}
/******
 *
 *  scatter
 *
 *
 ********/
TEST_CASE("Schur::PatchIfaceScatter<2> scatter local interfaces are copied")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      Vector<1> global_vector(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      auto local_vector = scatter.getNewLocalPatchIfaceVector();

      for (int i = 0; i < global_vector.getNumLocalPatches(); i++) {
        auto iface = iface_domain.getInterfaces()[i];
        auto local_data = global_vector.getComponentView(0, i);
        Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { local_data[coord] = iface->global_index + 1 + coord[0]; });
      }
      auto state = scatter.scatterStart(global_vector, *local_vector);
      for (auto iface : iface_domain.getInterfaces()) {
        auto local_data = local_vector->getComponentView(0, iface->local_index);
        Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(iface->global_index + 1 + coord[0])); });
      }
    }
  }
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatter")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      Vector<1> global_vector(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      auto local_vector = scatter.getNewLocalPatchIfaceVector();

      for (int i = 0; i < global_vector.getNumLocalPatches(); i++) {
        auto iface = iface_domain.getInterfaces()[i];
        auto local_data = global_vector.getComponentView(0, i);
        Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { local_data[coord] = iface->global_index + 1 + coord[0]; });
      }
      auto state = scatter.scatterStart(global_vector, *local_vector);
      scatter.scatterFinish(state, global_vector, *local_vector);
      for (auto piinfo : iface_domain.getPatchIfaceInfos()) {
        for (Side<2> s : Side<2>::getValues()) {
          if (piinfo->pinfo.hasNbr(s)) {
            auto iface_info = piinfo->getIfaceInfo(s);
            auto local_data = local_vector->getComponentView(0, iface_info->patch_local_index);
            Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(iface_info->global_index + 1 + coord[0])); });
          }
        }
      }
    }
  }
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatter twice")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      Vector<1> global_vector(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      auto local_vector = scatter.getNewLocalPatchIfaceVector();

      for (int i = 0; i < global_vector.getNumLocalPatches(); i++) {
        auto iface = iface_domain.getInterfaces()[i];
        auto local_data = global_vector.getComponentView(0, i);
        Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { local_data[coord] = iface->global_index + 1 + coord[0]; });
      }
      auto state = scatter.scatterStart(global_vector, *local_vector);
      scatter.scatterFinish(state, global_vector, *local_vector);
      state = scatter.scatterStart(global_vector, *local_vector);
      scatter.scatterFinish(state, global_vector, *local_vector);
      for (auto piinfo : iface_domain.getPatchIfaceInfos()) {
        for (Side<2> s : Side<2>::getValues()) {
          if (piinfo->pinfo.hasNbr(s)) {
            auto iface_info = piinfo->getIfaceInfo(s);
            auto local_data = local_vector->getComponentView(0, iface_info->patch_local_index);
            Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(iface_info->global_index + 1 + coord[0])); });
          }
        }
      }
    }
  }
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatter with local vector already filled")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 10 }) {

      DomainReader<2> domain_reader(mesh_file, { n, n }, 0);
      Domain<2> domain = domain_reader.getFinerDomain();
      Schur::InterfaceDomain<2> iface_domain(domain);

      Schur::PatchIfaceScatter<2> scatter(iface_domain);

      Vector<1> global_vector(domain.getCommunicator(), { n }, 1, iface_domain.getNumLocalInterfaces(), 0);
      auto local_vector = scatter.getNewLocalPatchIfaceVector();
      local_vector->setWithGhost(99);

      for (int i = 0; i < global_vector.getNumLocalPatches(); i++) {
        auto iface = iface_domain.getInterfaces()[i];
        auto local_data = global_vector.getComponentView(0, i);
        Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { local_data[coord] = iface->global_index + 1 + coord[0]; });
      }
      auto state = scatter.scatterStart(global_vector, *local_vector);
      scatter.scatterFinish(state, global_vector, *local_vector);
      for (auto piinfo : iface_domain.getPatchIfaceInfos()) {
        for (Side<2> s : Side<2>::getValues()) {
          if (piinfo->pinfo.hasNbr(s)) {
            auto iface_info = piinfo->getIfaceInfo(s);
            auto local_data = local_vector->getComponentView(0, iface_info->patch_local_index);
            Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(iface_info->global_index + 1 + coord[0])); });
          }
        }
      }
    }
  }
}
