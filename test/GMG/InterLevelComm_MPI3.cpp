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
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/InterLevelComm.h>

#include <doctest.h>

using namespace ThunderEgg;
using namespace std;

const string mesh_file = "mesh_inputs/2d_uniform_quad_mpi3.json";
TEST_CASE("3-processor InterLevelComm GetPatches on uniform quad")
{
  for (auto nx : { 2 }) {
    for (auto ny : { 2 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0) {
        CHECK_EQ(ilc.getPatchesWithGhostParent().size(), 0);
        CHECK_EQ(ilc.getPatchesWithLocalParent().size(), 2);

        map<int, set<int>> parents_to_children;
        for (auto pair : ilc.getPatchesWithLocalParent()) {
          parents_to_children[pair.first].insert(pair.second.get().id);
        }
        CHECK_UNARY(parents_to_children.count(0));
        CHECK_UNARY(parents_to_children[0].count(1));
        CHECK_UNARY(parents_to_children[0].count(2));
      } else if (rank == 1) {
        CHECK_EQ(ilc.getPatchesWithGhostParent().size(), 1);
        CHECK_EQ(ilc.getPatchesWithLocalParent().size(), 0);

        map<int, set<int>> parents_to_children;
        for (auto pair : ilc.getPatchesWithGhostParent()) {
          parents_to_children[pair.first].insert(pair.second.get().id);
        }
        CHECK_EQ(parents_to_children.count(0), 1);
        CHECK_EQ(parents_to_children[0].count(3), 1);
      } else {
        CHECK_EQ(ilc.getPatchesWithGhostParent().size(), 1);
        CHECK_EQ(ilc.getPatchesWithLocalParent().size(), 0);

        map<int, set<int>> parents_to_children;
        for (auto pair : ilc.getPatchesWithGhostParent()) {
          parents_to_children[pair.first].insert(pair.second.get().id);
        }
        CHECK_EQ(parents_to_children.count(0), 1);
        CHECK_EQ(parents_to_children[0].count(4), 1);
      }
    }
  }
}
TEST_CASE("3-processor getNewGhostVector on uniform quad")
{
  for (auto num_components : { 1, 2, 3 }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();
        Domain<2> d_coarse = domain_reader.getCoarserDomain();
        GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

        Vector<2> ghost_vec = ilc.getNewGhostVector(num_components);

        CHECK_EQ(ghost_vec.getNumComponents(), num_components);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
          CHECK_EQ(ghost_vec.getNumLocalPatches(), 0);
        } else {
          CHECK_EQ(ghost_vec.getNumLocalPatches(), 1);
        }
      }
    }
  }
}
TEST_CASE("3-processor sendGhostPatches on uniform quad")
{
  for (auto num_components : { 1, 2, 3 }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();
        Domain<2> d_coarse = domain_reader.getCoarserDomain();
        GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

        Vector<2> coarse_vec(d_coarse, num_components);

        Vector<2> ghost_vec = ilc.getNewGhostVector(num_components);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // fill vectors with rank+c+1
        for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
          PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
          Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) {
            for (int c = 0; c < num_components; c++) {
              local_view[coord] = rank + coord[2] + 1;
            }
          });
        }
        for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
          PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
          Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) {
            for (int c = 0; c < num_components; c++) {
              local_view[coord] = rank + coord[2] + 1;
            }
          });
        }

        ilc.sendGhostPatchesStart(coarse_vec, ghost_vec);
        ilc.sendGhostPatchesFinish(coarse_vec, ghost_vec);
        if (rank == 0) {
          // the coarse vec should be filled with 6+2*c
          PatchView<double, 2> local_view = coarse_vec.getPatchView(0);
          Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) {
            for (int c = 0; c < num_components; c++) {
              CHECK_EQ(local_view[coord], 1 + 2 + 3 + 3 * coord[2]);
            }
          });
        } else {
        }
      }
    }
  }
}
TEST_CASE("3-processor sendGhostPatches throws exception when start isn't called before finish on "
          "uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);

      CHECK_THROWS_AS(ilc.sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
    }
  }
}
TEST_CASE("3-processor sendGhostPatches throws exception when start and finish are called on "
          "different ghost vectors on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);
      Vector<2> ghost_vec_2 = ilc.getNewGhostVector(1);

      ilc.sendGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.sendGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
    }
  }
}
TEST_CASE("3-processor sendGhostPatches throws exception when start and finish are called on "
          "different vectors on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);
      Vector<2> coarse_vec_2(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);

      ilc.sendGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
    }
  }
}
TEST_CASE("3-processor sendGhostPatches throws exception when start and finish are called on "
          "different vectors and ghost vectors on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);
      Vector<2> coarse_vec_2(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);
      Vector<2> ghost_vec_2 = ilc.getNewGhostVector(1);

      ilc.sendGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.sendGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
    }
  }
}
TEST_CASE("3-processor sendGhostPatches throws exception when start is called twice on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);

      ilc.sendGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
    }
  }
}
TEST_CASE("3-processor getGhostPatches on uniform quad")
{
  for (auto num_components : { 1, 2, 3 }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();
        Domain<2> d_coarse = domain_reader.getCoarserDomain();
        GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

        Vector<2> coarse_vec(d_coarse, num_components);

        Vector<2> ghost_vec = ilc.getNewGhostVector(num_components);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // fill vectors with rank+c+1
        for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
          PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
          Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
        }
        for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
          PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
          Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
        }

        ilc.getGhostPatchesStart(coarse_vec, ghost_vec);
        ilc.getGhostPatchesFinish(coarse_vec, ghost_vec);
        if (rank == 0) {
        } else {
          // the coarse vec should be filled with 1+c
          PatchView<double, 2> local_view = ghost_vec.getPatchView(0);
          Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 1 + coord[2]); });
        }
      }
    }
  }
}
TEST_CASE("3-processor getGhostPatches throws exception when start isn't called before finish on "
          "uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);

      CHECK_THROWS_AS(ilc.getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
    }
  }
}
TEST_CASE("3-processor getGhostPatches throws exception when start and finish are called on "
          "different ghost vectors on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);
      Vector<2> ghost_vec_2 = ilc.getNewGhostVector(1);

      ilc.getGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.getGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
    }
  }
}
TEST_CASE("3-processor getGhostPatches throws exception when start and finish are called on "
          "different vectors on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);
      Vector<2> coarse_vec_2(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);

      ilc.getGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
    }
  }
}
TEST_CASE("3-processor getGhostPatches throws exception when start and finish are called on "
          "different vectors and ghost vectors on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);
      Vector<2> coarse_vec_2(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);
      Vector<2> ghost_vec_2 = ilc.getNewGhostVector(1);

      ilc.getGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.getGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
    }
  }
}
TEST_CASE("3-processor getGhostPatches throws exception when start is called twice on uniform quad")
{
  for (auto nx : { 2, 10 }) {
    for (auto ny : { 2, 10 }) {
      int num_ghost = 1;
      DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
      Domain<2> d_fine = domain_reader.getFinerDomain();
      Domain<2> d_coarse = domain_reader.getCoarserDomain();
      GMG::InterLevelComm<2> ilc(d_coarse, d_fine);

      Vector<2> coarse_vec(d_coarse, 1);

      Vector<2> ghost_vec = ilc.getNewGhostVector(1);

      ilc.getGhostPatchesStart(coarse_vec, ghost_vec);
      CHECK_THROWS_AS(ilc.getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
    }
  }
}
