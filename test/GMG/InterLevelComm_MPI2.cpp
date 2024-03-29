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

const string uniform = "mesh_inputs/2d_uniform_quad_mpi2.json";
const string mid_uniform = "mesh_inputs/2d_4x4_mid_on_1_mpi2.json";
#define MESHES uniform
#define MESHE_FILES uniform, mid_uniform

TEST_CASE("InterLevelComm Check number of local and ghost parents")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2 }) {
      for (auto ny : { 2 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          size_t num_ghost_parents = 0;
          size_t num_local_parents = 0;
          map<int, set<int>> my_local_parents_to_children;
          for (auto pinfo : d_fine.getPatchInfoVector()) {
            if (pinfo.parent_rank == rank) {
              num_local_parents++;
            } else {
              num_ghost_parents++;
            }
          }
          CHECK_EQ(ilc->getPatchesWithGhostParent().size(), num_ghost_parents);
          CHECK_EQ(ilc->getPatchesWithLocalParent().size(), num_local_parents);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm Check that parents have unique local indexes")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2 }) {
      for (auto ny : { 2 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);

          map<int, set<int>> local_index_id_map;
          for (auto pair : ilc->getPatchesWithLocalParent()) {
            local_index_id_map[pair.first].insert(pair.second.get().parent_id);
          }
          for (auto pair : local_index_id_map) {
            CHECK_EQ(pair.second.size(), 1);
          }
          map<int, set<int>> ghost_local_index_id_map;
          for (auto pair : ilc->getPatchesWithGhostParent()) {
            ghost_local_index_id_map[pair.first].insert(pair.second.get().parent_id);
          }
          for (auto pair : ghost_local_index_id_map) {
            CHECK_EQ(pair.second.size(), 1);
          }
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm Check that getPatches points to correct patches")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2 }) {
      for (auto ny : { 2 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);

          for (auto pair : ilc->getPatchesWithLocalParent()) {
            CHECK_EQ(&pair.second.get(), &ilc->getFinerDomain().getPatchInfoVector().at(pair.second.get().local_index));
          }
          for (auto pair : ilc->getPatchesWithGhostParent()) {
            CHECK_EQ(&pair.second.get(), &ilc->getFinerDomain().getPatchInfoVector().at(pair.second.get().local_index));
          }
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getNewGhostVector on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto num_components : { 1, 2, 3 }) {
      for (auto nx : { 2, 10 }) {
        for (auto ny : { 2, 10 }) {
          for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            Domain<2> d_coarse = domain_reader.getCoarserDomain();
            GMG::InterLevelComm<2>* ilc = nullptr;
            if (construction_type == "direct") {
              // direct construction
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            } else if (construction_type == "copy") {
              // copy construction
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(ilc_other);
            } else if (construction_type == "copy_assign") {
              // copy assignment
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
              *ilc = ilc_other;
            }

            Vector<2> ghost_vec = ilc->getNewGhostVector(num_components);

            CHECK_EQ(ghost_vec.getNumComponents(), num_components);
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0) {
              CHECK_EQ(ghost_vec.getNumLocalPatches(), 0);
            } else {
              CHECK_EQ(ghost_vec.getNumLocalPatches(), 1);
            }
            delete ilc;
          }
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto num_components : { 1, 2, 3 }) {
      for (auto nx : { 2, 10 }) {
        for (auto ny : { 2, 10 }) {
          for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            Domain<2> d_coarse = domain_reader.getCoarserDomain();
            GMG::InterLevelComm<2>* ilc = nullptr;
            if (construction_type == "direct") {
              // direct construction
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            } else if (construction_type == "copy") {
              // copy construction
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(ilc_other);
            } else if (construction_type == "copy_assign") {
              // copy assignment
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
              *ilc = ilc_other;
            }

            Vector<2> coarse_vec(d_coarse, num_components);

            Vector<2> ghost_vec = ilc->getNewGhostVector(num_components);

            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            // info

            // fill vectors with rank+c+1
            for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }
            for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }

            ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
            ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
            if (rank == 0) {
              // the coarse vec should be filled with 3+2*c
              PatchView<double, 2> local_view = coarse_vec.getPatchView(0);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 3 + 2 * coord[2]); });
            } else {
            }
            delete ilc;
          }
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches throws exception when start isn't called "
          "before finish on "
          "uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2 }) {
      for (auto ny : { 2 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches throws exception when start and finish are "
          "called on "
          "different ghost vectors on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2 }) {
      for (auto ny : { 2 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);
          Vector<2> ghost_vec_2 = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches throws exception when start and finish are "
          "called on "
          "different vectors on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);
          Vector<2> coarse_vec_2(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches throws exception when start and finish are "
          "called on "
          "different vectors and ghost vectors on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);
          Vector<2> coarse_vec_2(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);
          Vector<2> ghost_vec_2 = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("2-processor sendGhostPatches throws exception when start is called twice on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches throws exception when get start is called "
          "after send start "
          "on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches throws exception when sned start is called "
          "after get start "
          "on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto num_components : { 1, 2, 3 }) {
      for (auto nx : { 2, 10 }) {
        for (auto ny : { 2, 10 }) {
          for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            Domain<2> d_coarse = domain_reader.getCoarserDomain();
            GMG::InterLevelComm<2>* ilc = nullptr;
            if (construction_type == "direct") {
              // direct construction
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            } else if (construction_type == "copy") {
              // copy construction
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(ilc_other);
            } else if (construction_type == "copy_assign") {
              // copy assignment
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
              *ilc = ilc_other;
            }

            Vector<2> coarse_vec(d_coarse, num_components);

            Vector<2> ghost_vec = ilc->getNewGhostVector(1);

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

            ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
            ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
            if (rank == 0) {
            } else {
              // the coarse vec should be filled with 1+c
              PatchView<double, 2> local_view = ghost_vec.getPatchView(0);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 1 + coord[2]); });
            }
            delete ilc;
          }
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when start isn't called "
          "before finish on "
          "uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when start and finish are called on "
          "different ghost vectors on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);
          Vector<2> ghost_vec_2 = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when start and finish are called on "
          "different vectors on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);
          Vector<2> coarse_vec_2(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when start and finish are called on "
          "different vectors and ghost vectors on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);
          Vector<2> coarse_vec_2(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);
          Vector<2> ghost_vec_2 = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when start is called twice "
          "on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when send finish is called "
          "on get start")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when get finish is called "
          "on send start")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when send start is called "
          "after get start")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when get start is called "
          "after send start")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when send finish is called "
          "after get start "
          "and finish")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
          ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches throws exception when get finish is called "
          "after send start "
          "and finish")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();
          Domain<2> d_coarse = domain_reader.getCoarserDomain();
          GMG::InterLevelComm<2>* ilc = nullptr;
          if (construction_type == "direct") {
            // direct construction
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
          } else if (construction_type == "copy") {
            // copy construction
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(ilc_other);
          } else if (construction_type == "copy_assign") {
            // copy assignment
            GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
            ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            *ilc = ilc_other;
          }

          Vector<2> coarse_vec(d_coarse, 1);

          Vector<2> ghost_vec = ilc->getNewGhostVector(1);

          ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
          ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
          CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
          delete ilc;
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches called twice on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto num_components : { 1, 2, 3 }) {
      for (auto nx : { 2, 10 }) {
        for (auto ny : { 2, 10 }) {
          for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            Domain<2> d_coarse = domain_reader.getCoarserDomain();
            GMG::InterLevelComm<2>* ilc = nullptr;
            if (construction_type == "direct") {
              // direct construction
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            } else if (construction_type == "copy") {
              // copy construction
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(ilc_other);
            } else if (construction_type == "copy_assign") {
              // copy assignment
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
              *ilc = ilc_other;
            }

            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            for (int i = 0; i < 2; i++) {
              Vector<2> coarse_vec(d_coarse, num_components);

              Vector<2> ghost_vec = ilc->getNewGhostVector(num_components);

              // fill vectors with rank+c+1
              for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
                PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
                Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
              }
              for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
                PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
                Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
              }

              ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
              ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
              if (rank == 0) {
              } else {
                // the coarse vec should be filled with 1+c
                PatchView<double, 2> local_view = ghost_vec.getPatchView(0);
                Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 1 + coord[2]); });
              }
            }
            delete ilc;
          }
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches called twice on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto num_components : { 1, 2, 3 }) {
      for (auto nx : { 2, 10 }) {
        for (auto ny : { 2, 10 }) {
          for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            Domain<2> d_coarse = domain_reader.getCoarserDomain();
            GMG::InterLevelComm<2>* ilc = nullptr;
            if (construction_type == "direct") {
              // direct construction
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            } else if (construction_type == "copy") {
              // copy construction
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(ilc_other);
            } else if (construction_type == "copy_assign") {
              // copy assignment
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
              *ilc = ilc_other;
            }

            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            for (int i = 0; i < 2; i++) {
              Vector<2> coarse_vec(d_coarse, num_components);

              Vector<2> ghost_vec = ilc->getNewGhostVector(num_components);

              // fill vectors with rank+c+1
              for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
                PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
                Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
              }
              for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
                PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
                Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
              }

              ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
              ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
              if (rank == 0) {
                // the coarse vec should be filled with 3+2*c
                PatchView<double, 2> local_view = coarse_vec.getPatchView(0);
                Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 3 + 2 * coord[2]); });
              } else {
              }
            }
            delete ilc;
          }
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor sendGhostPatches then getGhostPaches called on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto num_components : { 1, 2, 3 }) {
      for (auto nx : { 2, 10 }) {
        for (auto ny : { 2, 10 }) {
          for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            Domain<2> d_coarse = domain_reader.getCoarserDomain();
            GMG::InterLevelComm<2>* ilc = nullptr;
            if (construction_type == "direct") {
              // direct construction
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            } else if (construction_type == "copy") {
              // copy construction
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(ilc_other);
            } else if (construction_type == "copy_assign") {
              // copy assignment
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
              *ilc = ilc_other;
            }

            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            Vector<2> coarse_vec(d_coarse, num_components);
            auto ghost_vec = ilc->getNewGhostVector(num_components);

            // fill vectors with rank+c+1
            for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }
            for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }

            ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
            ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
            if (rank == 0) {
              // the coarse vec should be filled with 3+2*c
              PatchView<double, 2> local_view = coarse_vec.getPatchView(0);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 3 + 2 * coord[2]); });
            } else {
            }

            // fill vectors with rank+c+1
            for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }
            for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }

            ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
            ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
            if (rank == 0) {
            } else {
              // the coarse vec should be filled with 1+c
              PatchView<double, 2> local_view = ghost_vec.getPatchView(0);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 1 + coord[2]); });
            }
            delete ilc;
          }
        }
      }
    }
  }
}
TEST_CASE("InterLevelComm 2-processor getGhostPatches then sendGhostPaches called on uniform quad")
{
  for (auto mesh_file : { MESHES }) {
    for (auto num_components : { 1, 2, 3 }) {
      for (auto nx : { 2, 10 }) {
        for (auto ny : { 2, 10 }) {
          for (std::string construction_type : { "direct", "copy", "copy_assign" }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            Domain<2> d_coarse = domain_reader.getCoarserDomain();
            GMG::InterLevelComm<2>* ilc = nullptr;
            if (construction_type == "direct") {
              // direct construction
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
            } else if (construction_type == "copy") {
              // copy construction
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(ilc_other);
            } else if (construction_type == "copy_assign") {
              // copy assignment
              GMG::InterLevelComm<2> ilc_other(d_coarse, d_fine);
              ilc = new GMG::InterLevelComm<2>(d_coarse, d_fine);
              *ilc = ilc_other;
            }

            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            Vector<2> coarse_vec(d_coarse, num_components);
            auto ghost_vec = ilc->getNewGhostVector(num_components);

            // fill vectors with rank+c+1
            for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }
            for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }

            ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
            ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
            if (rank == 0) {
            } else {
              // the coarse vec should be filled with 1+c
              PatchView<double, 2> local_view = ghost_vec.getPatchView(0);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 1 + coord[2]); });
            }

            // fill vectors with rank+c+1
            for (int i = 0; i < coarse_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = coarse_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }
            for (int i = 0; i < ghost_vec.getNumLocalPatches(); i++) {
              PatchView<double, 2> local_view = ghost_vec.getPatchView(i);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { local_view[coord] = rank + coord[2] + 1; });
            }

            ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
            ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
            if (rank == 0) {
              // the coarse vec should be filled with 3+2*c
              PatchView<double, 2> local_view = coarse_vec.getPatchView(0);
              Loop::OverAllIndexes<3>(local_view, [&](const std::array<int, 3>& coord) { CHECK_EQ(local_view[coord], 3 + 2 * coord[2]); });
            } else {
            }
            delete ilc;
          }
        }
      }
    }
  }
}
