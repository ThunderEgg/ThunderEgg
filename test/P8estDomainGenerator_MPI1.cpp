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

#include <doctest.h>

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
      p8est, false, [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return 1; }, nullptr);
    p8est_refine(
      p8est, false, [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return 1; }, nullptr);

    p8est_partition(p8est, true, nullptr);
  }
  ~FourTreeBSW()
  {
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);
  }
};
} // namespace
TEST_CASE("P8estDomainGenerator 4x4x4 hasCoarserDomain")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_2 = dg.getFinestDomain();
  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_1 = dg.getCoarserDomain();
  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_0 = dg.getCoarserDomain();
  CHECK_UNARY_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform Number of Patches")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CHECK_EQ(domain_2.getNumGlobalPatches(), 64);
  CHECK_EQ(domain_1.getNumGlobalPatches(), 8);
  CHECK_EQ(domain_0.getNumGlobalPatches(), 1);
}
TEST_CASE("P8estDomainGenerator 4x4x4 RefineLevel")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 2);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 1);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 0);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 rank")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 spacings")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (int nz : { 5, 10 }) {
        for (double scale_x : { 0.5, 1.0 }) {
          for (double scale_y : { 0.5, 1.0 }) {
            for (double scale_z : { 0.5, 1.0 }) {
              FourTreeBSW tree;

              P8estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
                x = scale_x * unit_x;
                y = scale_y * unit_y;
                z = scale_z * unit_z;
              };

              P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, bmf);

              Domain<3> domain_2 = dg.getFinestDomain();
              Domain<3> domain_1 = dg.getCoarserDomain();
              Domain<3> domain_0 = dg.getCoarserDomain();

              for (auto patch : domain_2.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.25 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.25 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.25 / nz));
              }

              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.5 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.5 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.5 / nz));
              }

              for (auto patch : domain_0.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 1.0 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 1.0 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 1.0 / nz));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 ns")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (int nz : { 5, 10 }) {
        FourTreeBSW tree;

        P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, Ident);

        Domain<3> domain_2 = dg.getFinestDomain();
        Domain<3> domain_1 = dg.getCoarserDomain();
        Domain<3> domain_0 = dg.getCoarserDomain();

        for (auto patch : domain_2.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }

        for (auto patch : domain_1.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }

        for (auto patch : domain_0.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 num_ghost_cells")
{
  for (int num_ghost_cells : { 0, 1, 2 }) {
    FourTreeBSW tree;

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, Ident);

    Domain<3> domain_2 = dg.getFinestDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();
    Domain<3> domain_0 = dg.getCoarserDomain();

    for (auto patch : domain_2.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }

    for (auto patch : domain_1.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }

    for (auto patch : domain_0.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform neighbor nfos")
{
  FourTreeBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CheckRootDomainNeighbors(domain_0);

  Check2x2x2DomainNeighbors(domain_1);

  Check4x4x4DomainNeighbors(domain_2);
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform child/parent ids and ranks")
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
      p8est, false, [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return 1; }, nullptr);
    p8est_refine(
      p8est, false, [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return 1; }, nullptr);
    p8est_refine(
      p8est, false, [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0; }, nullptr);

    p8est_partition(p8est, true, nullptr);
  }
  ~FourTreeRefineBSW()
  {
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);
  }
};
} // namespace
TEST_CASE("P8estDomainGenerator 4x4x4rbsw hasCoarserDomain")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_3 = dg.getFinestDomain();
  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_2 = dg.getCoarserDomain();
  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_1 = dg.getCoarserDomain();
  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_0 = dg.getCoarserDomain();
  CHECK_UNARY_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw Uniform Number of Patches")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CHECK_EQ(domain_3.getNumGlobalPatches(), 71);
  CHECK_EQ(domain_2.getNumGlobalPatches(), 64);
  CHECK_EQ(domain_1.getNumGlobalPatches(), 8);
  CHECK_EQ(domain_0.getNumGlobalPatches(), 1);
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw RefineLevel")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_3.getPatchInfoVector()) {
    if (patch.starts[0] < 0.24 && patch.starts[1] < 0.24 && patch.starts[2] < 0.24) {
      CHECK_EQ(patch.refine_level, 3);
    } else {
      CHECK_EQ(patch.refine_level, 2);
    }
  }

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 2);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 1);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 0);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw rank")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_3.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw spacings")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (int nz : { 5, 10 }) {
        for (double scale_x : { 0.5, 1.0 }) {
          for (double scale_y : { 0.5, 1.0 }) {
            for (double scale_z : { 0.5, 1.0 }) {
              FourTreeRefineBSW tree;

              P8estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
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
                if (patch.starts[0] < 0.24 * scale_x && patch.starts[1] < 0.24 * scale_y && patch.starts[2] < 0.24 * scale_z) {
                  CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.125 / nx));
                  CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.125 / ny));
                  CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.125 / nz));
                } else {
                  CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.25 / nx));
                  CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.25 / ny));
                  CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.25 / nz));
                }
              }

              for (auto patch : domain_2.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.25 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.25 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.25 / nz));
              }

              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.5 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.5 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.5 / nz));
              }

              for (auto patch : domain_0.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 1.0 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 1.0 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 1.0 / nz));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw ns")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (int nz : { 5, 10 }) {
        FourTreeRefineBSW tree;

        P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, Ident);

        Domain<3> domain_3 = dg.getFinestDomain();
        Domain<3> domain_2 = dg.getCoarserDomain();
        Domain<3> domain_1 = dg.getCoarserDomain();
        Domain<3> domain_0 = dg.getCoarserDomain();

        for (auto patch : domain_3.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }

        for (auto patch : domain_2.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }

        for (auto patch : domain_1.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }

        for (auto patch : domain_0.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw num_ghost_cells")
{
  for (int num_ghost_cells : { 0, 1, 2 }) {
    FourTreeRefineBSW tree;

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, Ident);

    Domain<3> domain_3 = dg.getFinestDomain();
    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();
    Domain<3> domain_0 = dg.getCoarserDomain();

    for (auto patch : domain_3.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }

    for (auto patch : domain_2.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }

    for (auto patch : domain_1.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }

    for (auto patch : domain_0.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw neighbor nfos")
{
  FourTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_3 = dg.getFinestDomain();
  Domain<3> domain_2 = dg.getCoarserDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CheckRootDomainNeighbors(domain_0);

  Check2x2x2DomainNeighbors(domain_1);

  Check4x4x4DomainNeighbors(domain_2);

  Check4x4x4RefinedBSWDomainNeighbors(domain_3);
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw child/parent ids and ranks")
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
      p8est, false, [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return 1; }, nullptr);
    p8est_refine(
      p8est, false, [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0; }, nullptr);

    p8est_partition(p8est, true, nullptr);
  }
  ~TwoTreeRefineBSW()
  {
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);
  }
};

} // namespace

TEST_CASE("P8estDomainGenerator 2x2x2rbsw Uniform hasCoarserDomain")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_2 = dg.getFinestDomain();
  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_1 = dg.getCoarserDomain();
  CHECK_UNARY(dg.hasCoarserDomain());
  Domain<3> domain_0 = dg.getCoarserDomain();
  CHECK_UNARY_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw Uniform Number of Patches")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CHECK_EQ(domain_2.getNumGlobalPatches(), 15);
  CHECK_EQ(domain_1.getNumGlobalPatches(), 8);
  CHECK_EQ(domain_0.getNumGlobalPatches(), 1);
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw RefineLevel")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    if (patch.starts[0] < 0.5 && patch.starts[1] < 0.5 && patch.starts[2] < 0.5) {
      CHECK_EQ(patch.refine_level, 2);
    } else {
      CHECK_EQ(patch.refine_level, 1);
    }
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 1);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK_EQ(patch.refine_level, 0);
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw rank")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  for (auto patch : domain_2.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }

  for (auto patch : domain_1.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }

  for (auto patch : domain_0.getPatchInfoVector()) {
    CHECK_EQ(patch.rank, tree.rank);
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw spacings")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (int nz : { 5, 10 }) {
        for (double scale_x : { 0.5, 1.0 }) {
          for (double scale_y : { 0.5, 1.0 }) {
            for (double scale_z : { 0.5, 1.0 }) {
              TwoTreeRefineBSW tree;

              P8estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
                x = scale_x * unit_x;
                y = scale_y * unit_y;
                z = scale_z * unit_z;
              };

              P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, bmf);

              Domain<3> domain_2 = dg.getFinestDomain();
              Domain<3> domain_1 = dg.getCoarserDomain();
              Domain<3> domain_0 = dg.getCoarserDomain();

              for (auto patch : domain_2.getPatchInfoVector()) {
                if (patch.starts[0] < 0.5 * scale_x && patch.starts[1] < 0.5 * scale_y && patch.starts[2] < 0.5 * scale_z) {
                  CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.25 / nx));
                  CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.25 / ny));
                  CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.25 / nz));
                } else {
                  CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.5 / nx));
                  CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.5 / ny));
                  CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.5 / nz));
                }
              }

              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.5 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.5 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.5 / nz));
              }

              for (auto patch : domain_0.getPatchInfoVector()) {
                CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 1.0 / nx));
                CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 1.0 / ny));
                CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 1.0 / nz));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw ns")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (int nz : { 5, 10 }) {
        TwoTreeRefineBSW tree;

        P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, Ident);

        Domain<3> domain_2 = dg.getFinestDomain();
        Domain<3> domain_1 = dg.getCoarserDomain();
        Domain<3> domain_0 = dg.getCoarserDomain();

        for (auto patch : domain_2.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }

        for (auto patch : domain_1.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }

        for (auto patch : domain_0.getPatchInfoVector()) {
          CHECK_EQ(patch.ns[0], nx);
          CHECK_EQ(patch.ns[1], ny);
          CHECK_EQ(patch.ns[2], nz);
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw num_ghost_cells")
{
  for (int num_ghost_cells : { 0, 1, 2 }) {
    TwoTreeRefineBSW tree;

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, Ident);

    Domain<3> domain_2 = dg.getFinestDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();
    Domain<3> domain_0 = dg.getCoarserDomain();

    for (auto patch : domain_2.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }

    for (auto patch : domain_1.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }

    for (auto patch : domain_0.getPatchInfoVector()) {
      CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
    }
  }
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw neighbor nfos")
{
  TwoTreeRefineBSW tree;

  P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, Ident);

  Domain<3> domain_2 = dg.getFinestDomain();
  Domain<3> domain_1 = dg.getCoarserDomain();
  Domain<3> domain_0 = dg.getCoarserDomain();

  CheckRootDomainNeighbors(domain_0);

  Check2x2x2DomainNeighbors(domain_1);

  Check2x2x2RefinedBSWDomainNeighbors(domain_2);
}
TEST_CASE("P8estDomainGenerator 2x2x2rbsw child/parent ids and ranks")
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
