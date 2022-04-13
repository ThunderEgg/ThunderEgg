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

#include "P8estDomainGenerator_SHARED.h"
#include <ThunderEgg/P8estDomainGenerator.h>

#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_mesh.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

namespace {
struct FourTreeBSW
{
  p8est_connectivity_t* conn;
  p8est_t* p8est;
  int rank;
  FourTreeBSW()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    conn = p8est_connectivity_new_unitcube();

    p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return 1;
      },
      nullptr);
    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return 1;
      },
      nullptr);

    p8est_partition(p8est, true, nullptr);
  }
  ~FourTreeBSW()
  {
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);
  }
};
} // namespace
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform hasCoarserDomain", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_2 = dg.getFinestDomain();
  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_1 = dg.getCoarserDomain();
  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_0 = dg.getCoarserDomain();
  CHECK_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform Number of Patches", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CHECK(domain_2.getNumGlobalPatches() == 64);
  CHECK(domain_1.getNumGlobalPatches() == 8);
  CHECK(domain_0.getNumGlobalPatches() == 1);
}
TEST_CASE("P8estDomainGenerator 4x4x4 RefineLevel", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.refine_level == 2);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.refine_level == 1);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.refine_level == 0);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 rank", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 spacings", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  int nx = GENERATE(5, 10);
  int ny = GENERATE(5, 10);
  int nz = GENERATE(5, 10);
  double scale_x = GENERATE(0.5, 1.0);
  double scale_y = GENERATE(0.5, 1.0);
  double scale_z = GENERATE(0.5, 1.0);

  P8estDomainGenerator::BlockMapFunc bmf =
    [&](
      int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
      x = scale_x * unit_x;
      y = scale_y * unit_y;
      z = scale_z * unit_z;
    };

  P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, bmf);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.25 / nz));
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.5 / nz));
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 1.0 / nz));
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 ns", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  int nx = GENERATE(5, 10);
  int ny = GENERATE(5, 10);
  int nz = GENERATE(5, 10);

  P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 num_ghost_cells", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  int num_ghost_cells = GENERATE(0, 1, 2);

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform neighbor nfos", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  INFO("Checking domain 0 neighbors");

  CheckRootDomainNeighbors(domain_0);

  INFO("Checking domain 1 neighbors");

  Check2x2x2DomainNeighbors(domain_1);

  INFO("Checking domain 2 neighbors");

  Check4x4x4DomainNeighbors(domain_2);
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform child/parent ids and ranks", "[P8estDomainGenerator]")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CheckChildIdsAndRanksNull(domain_2);

  CheckParentAndChildIdsAndRanks(domain_1, 1, domain_2, 2);

  CheckParentAndChildIdsAndRanks(domain_0, 0, domain_1, 1);

  CheckParentIdsAndRanksNull(domain_0);

  PatchVector domain_2_pvector(domain_2, 2);
  PatchVector domain_1_pvector(domain_1, 1);
  PatchVector domain_0_pvector(domain_0, 0);
}
namespace {
struct FourTreeRefineBSW
{
  p8est_connectivity_t* conn;
  p8est_t* p8est;
  int rank;
  FourTreeRefineBSW()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    conn = p8est_connectivity_new_unitcube();

    p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return 1;
      },
      nullptr);
    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return 1;
      },
      nullptr);
    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0;
      },
      nullptr);

    p8est_partition(p8est, true, nullptr);
  }
  ~FourTreeRefineBSW()
  {
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);
  }
};
} // namespace
TEST_CASE("P8estDomainGenerator 4x4x4rbsw Uniform hasCoarserDomain", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_3 = dg.getFinestDomain();
  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_2 = dg.getCoarserDomain();
  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_1 = dg.getCoarserDomain();
  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_0 = dg.getCoarserDomain();
  CHECK_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw Uniform Number of Patches", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CHECK(domain_3.getNumGlobalPatches() == 71);
  CHECK(domain_2.getNumGlobalPatches() == 64);
  CHECK(domain_1.getNumGlobalPatches() == 8);
  CHECK(domain_0.getNumGlobalPatches() == 1);
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw RefineLevel", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_3.getPatchInfoVector()) {
    if (patch.starts[0] < 0.24 && patch.starts[1] < 0.24 && patch.starts[2] < 0.24) {
      CHECK(patch.refine_level == 3);
    } else {
      CHECK(patch.refine_level == 2);
    }
  }

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.refine_level == 2);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.refine_level == 1);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.refine_level == 0);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw rank", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_3.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw spacings", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  int nx = GENERATE(5, 10);
  int ny = GENERATE(5, 10);
  int nz = GENERATE(5, 10);
  double scale_x = GENERATE(0.5, 1.0);
  double scale_y = GENERATE(0.5, 1.0);
  double scale_z = GENERATE(0.5, 1.0);

  P8estDomainGenerator::BlockMapFunc bmf =
    [&](
      int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
      x = scale_x * unit_x;
      y = scale_y * unit_y;
      z = scale_z * unit_z;
    };

  P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, bmf);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_3.getPatchInfoVector()) {
    if (patch.starts[0] < 0.24 * scale_x && patch.starts[1] < 0.24 * scale_y &&
        patch.starts[2] < 0.24 * scale_z) {
      CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.125 / nx));
      CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.125 / ny));
      CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.125 / nz));
    } else {
      CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
      CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
      CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.25 / nz));
    }
  }

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.25 / nz));
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.5 / nz));
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 1.0 / nz));
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw ns", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  int nx = GENERATE(5, 10);
  int ny = GENERATE(5, 10);
  int nz = GENERATE(5, 10);

  P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_3.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw num_ghost_cells", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  int num_ghost_cells = GENERATE(0, 1, 2);

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_3.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw neighbor nfos", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  INFO("Checking domain 0 neighbors");

  CheckRootDomainNeighbors(domain_0);

  INFO("Checking domain 1 neighbors");

  Check2x2x2DomainNeighbors(domain_1);

  INFO("Checking domain 2 neighbors");

  Check4x4x4DomainNeighbors(domain_2);

  INFO("Checking domain 3 neighbors");

  Check4x4x4RefinedBSWDomainNeighbors(domain_3);
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw child/parent ids and ranks", "[P8estDomainGenerator]")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CheckChildIdsAndRanksNull(domain_3);

  CheckParentAndChildIdsAndRanksRefined(domain_2, 2, domain_3, 3);

  CheckParentAndChildIdsAndRanks(domain_1, 1, domain_2, 2);

  CheckParentAndChildIdsAndRanks(domain_0, 0, domain_1, 1);

  CheckParentIdsAndRanksNull(domain_0);
}
namespace {
struct TwoTreeRefineBSW
{
  p8est_connectivity_t* conn;
  p8est_t* p8est;
  int rank;
  TwoTreeRefineBSW()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    conn = p8est_connectivity_new_unitcube();

    p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return 1;
      },
      nullptr);
    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0;
      },
      nullptr);

    p8est_partition(p8est, true, nullptr);
  }
  ~TwoTreeRefineBSW()
  {
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);
  }
};

} // namespace
TEST_CASE("P8estDomainGenerator 2x2x2rbsw Uniform hasCoarserDomain", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_2 = dg.getFinestDomain();
  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_1 = dg.getCoarserDomain();
  CHECK(dg.hasCoarserDomain());
  Domain<3> domain_0 = dg.getCoarserDomain();
  CHECK_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw Uniform Number of Patches", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CHECK(domain_2.getNumGlobalPatches() == 15);
  CHECK(domain_1.getNumGlobalPatches() == 8);
  CHECK(domain_0.getNumGlobalPatches() == 1);
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw RefineLevel", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    if (patch.starts[0] < 0.5 && patch.starts[1] < 0.5 && patch.starts[2] < 0.5) {
      CHECK(patch.refine_level == 2);
    } else {
      CHECK(patch.refine_level == 1);
    }
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.refine_level == 1);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.refine_level == 0);
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw rank", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.rank == tree.rank);
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw spacings", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  int nx = GENERATE(5, 10);
  int ny = GENERATE(5, 10);
  int nz = GENERATE(5, 10);
  double scale_x = GENERATE(0.5, 1.0);
  double scale_y = GENERATE(0.5, 1.0);
  double scale_z = GENERATE(0.5, 1.0);

  P8estDomainGenerator::BlockMapFunc bmf =
    [&](
      int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
      x = scale_x * unit_x;
      y = scale_y * unit_y;
      z = scale_z * unit_z;
    };

  P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, bmf);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    if (patch.starts[0] < 0.5 * scale_x && patch.starts[1] < 0.5 * scale_y &&
        patch.starts[2] < 0.5 * scale_z) {
      CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
      CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
      CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.25 / nz));
    } else {
      CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
      CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
      CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.5 / nz));
    }
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.5 / nz));
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
    CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
    CHECK(patch.spacings[2] == Catch::Approx(scale_z * 1.0 / nz));
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw ns", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  int nx = GENERATE(5, 10);
  int ny = GENERATE(5, 10);
  int nz = GENERATE(5, 10);

  P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.ns[0] == nx);
    CHECK(patch.ns[1] == ny);
    CHECK(patch.ns[2] == nz);
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw num_ghost_cells", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  int num_ghost_cells = GENERATE(0, 1, 2);

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK(patch.num_ghost_cells == num_ghost_cells);
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw neighbor nfos", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  INFO("Checking domain 0 neighbors");

  CheckRootDomainNeighbors(domain_0);

  INFO("Checking domain 1 neighbors");

  Check2x2x2DomainNeighbors(domain_1);

  INFO("Checking domain 2 neighbors");

  Check2x2x2RefinedBSWDomainNeighbors(domain_2);
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw child/parent ids and ranks", "[P8estDomainGenerator]")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CheckChildIdsAndRanksNull(domain_2);

  CheckParentAndChildIdsAndRanksRefined(domain_1, 1, domain_2, 2);

  CheckParentAndChildIdsAndRanks(domain_0, 0, domain_1, 1);

  CheckParentIdsAndRanksNull(domain_0);

  PatchVector domain_2_pvector(domain_2, 2);
  PatchVector domain_1_pvector(domain_1, 1);
  PatchVector domain_0_pvector(domain_0, 0);
}
